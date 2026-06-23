library(bsseq)
library(HDF5Array)
library(GenomicRanges)
library(limma)
library(TxDb.Mmulatta.UCSC.rheMac10.refGene)
library(org.Mmu.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

set.seed(42)

AGE_SPLIT    <- 10     # inflection point for Cluster 3+4 stratification (years)
UPSTREAM_KB  <- 2000   # bp upstream of TSS to include

# ============================================================
# 1. LOAD DATA
# ============================================================

cat("Loading BSseq data...\n")
bs <- loadHDF5SummarizedExperiment(
  "/home/tander10/sternerlab/brain_methylation/results/bs_filtered_86samples_hdf5"
)

cat("Loading metadata...\n")
metadata <- readRDS(
  "/home/tander10/sternerlab/brain_methylation/results/metadata_filtered_86samples_with_SV1-3.rds"
)
colnames(metadata)[colnames(metadata) == "Age"] <- "age"
colnames(metadata)[colnames(metadata) == "Sex"] <- "sex"
metadata$sex    <- as.factor(metadata$sex)
metadata$sampid <- trimws(as.character(metadata$sampid))

cat("Loading RNAseq residualized counts...\n")
rnaseq_resid <- read.csv(
  "/home/tander10/sternerlab/brain_methylation/RNAseq_residualizedcounts.csv",
  row.names = 1, check.names = FALSE
)
colnames(rnaseq_resid) <- trimws(as.character(colnames(rnaseq_resid)))

cat("Loading gene cluster assignments (RNAseq_genes.csv)...\n")
gene_clusters <- read.csv("RNAseq_genes.csv", stringsAsFactors = FALSE)

# Normalize column names
colnames(gene_clusters) <- trimws(tolower(colnames(gene_clusters)))
if (!"gene"    %in% colnames(gene_clusters)) colnames(gene_clusters)[1] <- "gene"
if (!"cluster" %in% colnames(gene_clusters)) colnames(gene_clusters)[2] <- "cluster"
gene_clusters$gene <- trimws(as.character(gene_clusters$gene))

cat("  Total genes in cluster file:", nrow(gene_clusters), "\n")
cat("  Cluster breakdown:\n")
print(table(gene_clusters$cluster))
cat("\n")

# Clusters 1+2 = linear trajectories; 3+4 = nonlinear (inflection ~10 yr)
LINEAR_CLUSTERS    <- c(1, 2)
NONLINEAR_CLUSTERS <- c(3, 4)

linear_genes     <- gene_clusters$gene[gene_clusters$cluster %in% LINEAR_CLUSTERS]
nonlinear_genes  <- gene_clusters$gene[gene_clusters$cluster %in% NONLINEAR_CLUSTERS]
gene_cluster_map <- setNames(gene_clusters$cluster, gene_clusters$gene)

cat("  Linear trajectory genes (Clusters 1+2):    ", length(linear_genes), "\n")
cat("  Nonlinear trajectory genes (Clusters 3+4): ", length(nonlinear_genes), "\n\n")

# ============================================================
# 2. MATCH SAMPLES ACROSS DNAm AND RNAseq DATASETS
# ============================================================

shared_samples <- intersect(metadata$sampid, colnames(rnaseq_resid))
cat("Shared samples (DNAm ∩ RNAseq):", length(shared_samples), "\n")

if (length(shared_samples) < 10) {
  stop("Too few shared samples — check that sampid matches RNAseq column names.")
}

metadata_shared <- metadata[match(shared_samples, metadata$sampid), ]
stopifnot(all(metadata_shared$sampid == shared_samples))

# Age-stratified subsets for nonlinear analysis
early_mask    <- metadata_shared$age <  AGE_SPLIT
late_mask     <- metadata_shared$age >= AGE_SPLIT
early_samples <- shared_samples[early_mask]
late_samples  <- shared_samples[late_mask]

cat(sprintf("  Age <  %d years (early):  %d samples | age range: %.1f–%.1f yr\n",
            AGE_SPLIT, sum(early_mask),
            min(metadata_shared$age[early_mask]),
            max(metadata_shared$age[early_mask])))
cat(sprintf("  Age >= %d years (late):   %d samples | age range: %.1f–%.1f yr\n\n",
            AGE_SPLIT, sum(late_mask),
            min(metadata_shared$age[late_mask]),
            max(metadata_shared$age[late_mask])))

if (sum(early_mask) < 5 || sum(late_mask) < 5) {
  warning("One age group has < 5 samples — models will have very low power. Adjust AGE_SPLIT.")
}

# ============================================================
# 3. EXTRACT CpG METHYLATION AND REGRESS OUT COVARIATES
# ============================================================

cat("Extracting raw methylation matrix...\n")
proportion_meth <- as.matrix(getMeth(bs, type = "raw"))

cpg_chr <- as.character(seqnames(rowRanges(bs)))
cpg_pos <- start(rowRanges(bs))
rownames(proportion_meth) <- paste0(cpg_chr, ":", cpg_pos)
colnames(proportion_meth) <- trimws(as.character(metadata$sampid))

proportion_meth <- proportion_meth[, shared_samples, drop = FALSE]
proportion_meth <- proportion_meth[complete.cases(proportion_meth), ]
cat("  CpGs retained after NA removal:", nrow(proportion_meth), "\n")

cat("Regressing out sex + SV1-SV3 from methylation (full dataset)...\n")
design_meth   <- model.matrix(~ sex + SV1 + SV2 + SV3, data = metadata_shared)
residual_meth <- removeBatchEffect(proportion_meth, design = design_meth)

# Build CpG GRanges from row names
retained_split <- strsplit(rownames(residual_meth), ":")
cpg_gr <- GRanges(
  seqnames = sapply(retained_split, `[`, 1),
  ranges   = IRanges(start = as.integer(sapply(retained_split, `[`, 2)), width = 1)
)

# ============================================================
# 4. DEFINE GENE BODY REGIONS: [TSS − 2kb] → [TES]
# ============================================================

cat("Building gene body regions (TSS −", UPSTREAM_KB, "bp → TES)...\n")
txdb     <- TxDb.Mmulatta.UCSC.rheMac10.refGene
genes_gr <- genes(txdb)

# Harmonise chromosome naming
cpg_style  <- seqlevelsStyle(cpg_gr)[1]
txdb_style <- seqlevelsStyle(genes_gr)[1]
if (cpg_style != txdb_style) {
  seqlevelsStyle(cpg_gr) <- txdb_style
  cat("  Chromosome styles harmonised:", cpg_style, "→", txdb_style, "\n")
}

# All cluster genes present in the expression matrix
all_cluster_genes <- c(linear_genes, nonlinear_genes)
all_input_genes   <- intersect(all_cluster_genes, rownames(rnaseq_resid))
cat("  Cluster genes present in expression matrix:", length(all_input_genes),
    "/", length(all_cluster_genes), "\n")

if (length(all_input_genes) == 0) {
  cat("\n--- DIAGNOSTIC ---\n")
  cat("First 10 gene names from RNAseq_genes.csv:\n")
  print(head(all_cluster_genes, 10))
  cat("First 10 row names from expression matrix:\n")
  print(head(rownames(rnaseq_resid), 10))
  stop(paste(
    "No genes matched. Check that gene symbols in RNAseq_genes.csv",
    "match row names in the residualized expression matrix."
  ))
}

# Map gene symbols → Entrez IDs
all_entrez <- mapIds(
  org.Mmu.eg.db,
  keys      = all_input_genes,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
cat("  Symbols mapped to Entrez IDs:", sum(!is.na(all_entrez)),
    "/", length(all_input_genes), "\n")

all_entrez    <- all_entrez[!is.na(all_entrez)]
entrez_to_sym <- setNames(names(all_entrez), all_entrez)

# Subset TxDb to mapped genes
input_genes_gr         <- genes_gr[names(genes_gr) %in% names(entrez_to_sym)]
input_genes_gr$symbol  <- entrez_to_sym[names(input_genes_gr)]
input_genes_gr         <- input_genes_gr[!duplicated(input_genes_gr$symbol)]
cat("  Genes with Mmul10 coordinates:", length(input_genes_gr), "\n")

if (length(input_genes_gr) == 0) {
  stop("No genes could be mapped to genomic coordinates. Check gene symbol format.")
}

# --- Build gene body + upstream flank (strand-aware) ---
upstream_flank              <- trim(flank(input_genes_gr, width = UPSTREAM_KB, start = TRUE))
genebody_ext_gr             <- punion(input_genes_gr, upstream_flank, fill.gap = TRUE)
genebody_ext_gr$gene_symbol <- input_genes_gr$symbol

cat(sprintf("  Gene body regions defined: %d genes\n", length(genebody_ext_gr)))
cat(sprintf("  Median region width: %d bp (vs median gene body: %d bp)\n",
            median(width(genebody_ext_gr)),
            median(width(input_genes_gr))))

# ============================================================
# 5. COMPUTE MEAN GENE BODY METHYLATION PER GENE (ALL SAMPLES)
# ============================================================

cat("Computing mean gene body methylation per gene...\n")
hits    <- findOverlaps(genebody_ext_gr, cpg_gr)
n_cpgs  <- countOverlaps(genebody_ext_gr, cpg_gr)
has_cov <- n_cpgs > 0
cat("  Genes with ≥1 CpG in gene body region:", sum(has_cov),
    "/", length(genebody_ext_gr), "\n")

genes_covered <- genebody_ext_gr$gene_symbol[has_cov]
genebody_meth <- matrix(
  NA,
  nrow     = sum(has_cov),
  ncol     = length(shared_samples),
  dimnames = list(genes_covered, shared_samples)
)

for (i in which(has_cov)) {
  gene <- genebody_ext_gr$gene_symbol[i]
  idx  <- subjectHits(hits[queryHits(hits) == i])
  genebody_meth[gene, ] <- if (length(idx) == 1) {
    residual_meth[idx, ]
  } else {
    colMeans(residual_meth[idx, ], na.rm = TRUE)
  }
}

# Remove genes with any remaining NAs
complete_genes <- rowSums(is.na(genebody_meth)) == 0
genebody_meth  <- genebody_meth[complete_genes, ]
n_cpgs_vec     <- n_cpgs[has_cov][complete_genes]
names(n_cpgs_vec) <- rownames(genebody_meth)
cat("  Genes with complete gene body methylation:", nrow(genebody_meth), "\n\n")

# ============================================================
# 6. SUBSET MATRICES INTO LINEAR AND NONLINEAR GENE SETS
# ============================================================

linear_tested    <- intersect(rownames(genebody_meth),
                              intersect(rownames(rnaseq_resid), linear_genes))
nonlinear_tested <- intersect(rownames(genebody_meth),
                              intersect(rownames(rnaseq_resid), nonlinear_genes))

cat("Final genes entering models:\n")
cat("  Linear (Clusters 1+2):    ", length(linear_tested), "\n")
cat("  Nonlinear (Clusters 3+4): ", length(nonlinear_tested), "\n\n")

gb_linear      <- genebody_meth[linear_tested,    , drop = FALSE]
expr_linear    <- as.matrix(rnaseq_resid[linear_tested,    shared_samples, drop = FALSE])

gb_nonlinear   <- genebody_meth[nonlinear_tested, , drop = FALSE]
expr_nonlinear <- as.matrix(rnaseq_resid[nonlinear_tested, shared_samples, drop = FALSE])

# ============================================================
# HELPER — fit lm(expr ~ genebody_meth) across genes and samples
# ============================================================
run_genebody_lm <- function(gb_mat, expr_mat, samples, n_cpgs_vec, cluster_map) {
  genes <- intersect(rownames(gb_mat), rownames(expr_mat))

  res <- data.frame(
    gene       = genes,
    cluster    = as.character(cluster_map[genes]),
    n_cpgs     = as.integer(n_cpgs_vec[genes]),
    n_samples  = length(samples),
    beta_meth  = NA_real_,
    se_meth    = NA_real_,
    t_meth     = NA_real_,
    pval_meth  = NA_real_,
    r_squared  = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes)) {
    gene <- genes[i]
    meth <- as.numeric(gb_mat[gene,   samples])
    expr <- as.numeric(expr_mat[gene, samples])
    if (anyNA(meth) || anyNA(expr)) next

    tryCatch({
      fit <- lm(expr ~ meth)
      ct  <- summary(fit)$coefficients
      res$beta_meth[i] <- ct["meth", "Estimate"]
      res$se_meth[i]   <- ct["meth", "Std. Error"]
      res$t_meth[i]    <- ct["meth", "t value"]
      res$pval_meth[i] <- ct["meth", "Pr(>|t|)"]
      res$r_squared[i] <- summary(fit)$r.squared
    }, error = function(e) NULL)
  }

  res$fdr       <- p.adjust(res$pval_meth, method = "BH")
  res$direction <- ifelse(res$beta_meth < 0, "repressive", "activating")
  res$sig_05    <- !is.na(res$fdr) & res$fdr < 0.05
  res$sig_10    <- !is.na(res$fdr) & res$fdr < 0.10
  res
}

# ============================================================
# 7. ANALYSIS A — LINEAR GENES (Clusters 1+2), ALL SAMPLES
# ============================================================

cat("=== ANALYSIS A: Linear genes (Clusters 1+2) — all samples ===\n")
results_linear <- run_genebody_lm(
  gb_mat      = gb_linear,
  expr_mat    = expr_linear,
  samples     = shared_samples,
  n_cpgs_vec  = n_cpgs_vec,
  cluster_map = gene_cluster_map
)
cat("  Models fit:", nrow(results_linear), "\n")

# ============================================================
# 8. ANALYSIS B — NONLINEAR GENES (Clusters 3+4), STRATIFIED
# ============================================================

cat("\n=== ANALYSIS B: Nonlinear genes (Clusters 3+4) — age stratified ===\n")
cat(sprintf("  Inflection point: %d years\n", AGE_SPLIT))

cat(sprintf("\n  B-Early: age < %d (%d samples)...\n", AGE_SPLIT, sum(early_mask)))
results_early <- run_genebody_lm(
  gb_mat      = gb_nonlinear,
  expr_mat    = expr_nonlinear,
  samples     = early_samples,
  n_cpgs_vec  = n_cpgs_vec,
  cluster_map = gene_cluster_map
)

cat(sprintf("  B-Late: age >= %d (%d samples)...\n", AGE_SPLIT, sum(late_mask)))
results_late <- run_genebody_lm(
  gb_mat      = gb_nonlinear,
  expr_mat    = expr_nonlinear,
  samples     = late_samples,
  n_cpgs_vec  = n_cpgs_vec,
  cluster_map = gene_cluster_map
)

# ============================================================
# 9. SUMMARY STATISTICS
# ============================================================

print_summary <- function(results, label) {
  cat("\n---------------------------------------------------------------\n")
  cat(label, "\n")
  cat("---------------------------------------------------------------\n")
  cat("  Genes tested:", nrow(results),
      "| Samples:", unique(results$n_samples), "\n\n")

  for (fdr_cut in c(0.05, 0.10)) {
    sig_col <- if (fdr_cut == 0.05) "sig_05" else "sig_10"
    sig_res <- results[!is.na(results[[sig_col]]) & results[[sig_col]], ]
    cat(sprintf("  FDR < %.2f: %d genes", fdr_cut, nrow(sig_res)))
    if (nrow(sig_res) > 0) {
      cat(sprintf("  [repressive: %d (%.0f%%)  activating: %d (%.0f%%)]",
          sum(sig_res$direction == "repressive"),
          100 * mean(sig_res$direction == "repressive"),
          sum(sig_res$direction == "activating"),
          100 * mean(sig_res$direction == "activating")))
    }
    cat("\n")
  }

  n_all <- sum(!is.na(results$direction))
  n_rep <- sum(results$direction == "repressive", na.rm = TRUE)
  cat(sprintf("\n  Overall (%d genes): %d repressive (%.1f%%) | %d activating (%.1f%%)\n",
              n_all, n_rep, 100 * n_rep / n_all,
              n_all - n_rep, 100 * (n_all - n_rep) / n_all))
  binom_p <- binom.test(n_rep, n_all, p = 0.5, alternative = "greater")$p.value
  cat(sprintf("  Binomial test (repressive > 50%%): p = %.4f\n", binom_p))

  top5 <- results[!is.na(results$fdr), ]
  top5 <- top5[order(top5$fdr), ]
  top5 <- head(top5, 5)
  if (nrow(top5) > 0) {
    cat("\n  Top 5 by FDR:\n")
    print(top5[, c("gene", "cluster", "n_cpgs", "beta_meth",
                   "r_squared", "pval_meth", "fdr", "direction")],
          row.names = FALSE)
  }
}

cat("\n===============================================================\n")
cat("GENE BODY METHYLATION → GENE EXPRESSION — STRATIFIED RESULTS\n")
cat("Region: TSS −", UPSTREAM_KB, "bp → TES (full gene body + upstream flank)\n")
cat("===============================================================\n")

print_summary(results_linear,
  "ANALYSIS A — Linear trajectories (Clusters 1+2), all samples")
print_summary(results_early,
  sprintf("ANALYSIS B-EARLY — Nonlinear (Clusters 3+4), age < %d years", AGE_SPLIT))
print_summary(results_late,
  sprintf("ANALYSIS B-LATE  — Nonlinear (Clusters 3+4), age >= %d years", AGE_SPLIT))

# ============================================================
# 10. EARLY vs LATE COMPARISON (nonlinear genes)
# ============================================================

cat("\n===============================================================\n")
cat("EARLY vs LATE COMPARISON — Nonlinear Trajectory Genes\n")
cat("===============================================================\n")

comparison <- results_early %>%
  dplyr::select(gene, cluster, n_cpgs,
                beta_early = beta_meth, se_early  = se_meth,
                pval_early = pval_meth, fdr_early  = fdr,
                sig05_early = sig_05,  dir_early  = direction,
                r2_early    = r_squared) %>%
  dplyr::left_join(
    results_late %>%
      dplyr::select(gene,
                    beta_late  = beta_meth, se_late   = se_meth,
                    pval_late  = pval_meth, fdr_late   = fdr,
                    sig05_late = sig_05,   dir_late   = direction,
                    r2_late    = r_squared),
    by = "gene"
  ) %>%
  mutate(
    period_sig = case_when(
       sig05_early &  sig05_late ~ "Both periods",
       sig05_early & !sig05_late ~ "Early only",
      !sig05_early &  sig05_late ~ "Late only",
      TRUE                       ~ "Neither"
    ),
    direction_switch = !is.na(dir_early) & !is.na(dir_late) & (dir_early != dir_late),
    delta_beta       = beta_late - beta_early
  )

cat("\nSignificance pattern across life stages (FDR < 0.05):\n")
print(comparison %>% dplyr::count(period_sig, name = "n_genes"))

early_only   <- comparison %>% dplyr::filter(period_sig == "Early only") %>% dplyr::arrange(fdr_early)
late_only    <- comparison %>% dplyr::filter(period_sig == "Late only")  %>% dplyr::arrange(fdr_late)
both_sig     <- comparison %>%
  dplyr::filter(period_sig == "Both periods") %>%
  dplyr::mutate(avg_fdr = (fdr_early + fdr_late) / 2) %>%
  dplyr::arrange(avg_fdr)

cat("\n--- Genes significant in EARLY life only ---\n")
if (nrow(early_only) > 0) {
  print(early_only %>% dplyr::select(gene, cluster, beta_early, fdr_early,
                                      beta_late, fdr_late, direction_switch),
        row.names = FALSE)
} else { cat("  None at FDR < 0.05\n") }

cat("\n--- Genes significant in LATE life only ---\n")
if (nrow(late_only) > 0) {
  print(late_only %>% dplyr::select(gene, cluster, beta_early, fdr_early,
                                     beta_late, fdr_late, direction_switch),
        row.names = FALSE)
} else { cat("  None at FDR < 0.05\n") }

cat("\n--- Genes significant in BOTH periods ---\n")
if (nrow(both_sig) > 0) {
  print(both_sig %>% dplyr::select(gene, cluster, beta_early, fdr_early,
                                    beta_late, fdr_late, direction_switch, delta_beta),
        row.names = FALSE)
} else { cat("  None at FDR < 0.05\n") }

cat("\n--- Genes that SWITCH direction between periods ---\n")
switches <- comparison %>%
  dplyr::filter(direction_switch, period_sig != "Neither") %>%
  dplyr::arrange(fdr_early)
if (nrow(switches) > 0) {
  print(switches %>% dplyr::select(gene, cluster, dir_early, beta_early, fdr_early,
                                    dir_late, beta_late, fdr_late),
        row.names = FALSE)
} else { cat("  No direction switches among significant genes\n") }

# ============================================================
# 11. SAVE ALL RESULTS
# ============================================================

out_dir <- "/home/tander10/sternerlab/brain_methylation/genebody_GE_stratified"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  results_linear[order(results_linear$fdr, na.last = TRUE), ],
  file.path(out_dir, "AnalysisA_Linear_GenebodyMeth_Results.csv"),
  row.names = FALSE
)
write.csv(
  results_early[order(results_early$fdr, na.last = TRUE), ],
  file.path(out_dir, "AnalysisB_Nonlinear_Early_GenebodyMeth_Results.csv"),
  row.names = FALSE
)
write.csv(
  results_late[order(results_late$fdr, na.last = TRUE), ],
  file.path(out_dir, "AnalysisB_Nonlinear_Late_GenebodyMeth_Results.csv"),
  row.names = FALSE
)
write.csv(
  comparison[order(comparison$fdr_early, na.last = TRUE), ],
  file.path(out_dir, "AnalysisB_Nonlinear_EarlyLate_Comparison.csv"),
  row.names = FALSE
)
saveRDS(genebody_meth,
        file.path(out_dir, "genebody_meth_matrix_all_genes.rds"))
