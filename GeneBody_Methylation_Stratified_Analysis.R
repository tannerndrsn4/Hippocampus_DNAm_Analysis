# ============================================================
# GENE BODY METHYLATION → GENE EXPRESSION — STRATIFIED ANALYSIS
#
# Companion to Promoter_Methylation_Stratified_Analysis.R.
# Identical design but the CpG region covers the full gene body
# plus 2kb upstream of the TSS, rather than the promoter window
# (TSS ± 2kb).
#
#   Region per gene: [TSS − 2kb] → [TES]  (strand-aware)
#   + strand: genomic start − 2000  →  genomic end
#   − strand: genomic start          →  genomic end + 2000
#
# This tests whether methylation anywhere from just upstream of
# the promoter through the full gene body predicts expression.
#
# Three analyses:
#   Analysis A  — Linear trajectory genes (Clusters 1+2), all samples.
#   Analysis B-Early — Nonlinear genes (Clusters 3+4), age < 10 yr.
#   Analysis B-Late  — Nonlinear genes (Clusters 3+4), age >= 10 yr.
#
# Model: resid_expression ~ resid_genebody_meth
#   Both matrices pre-residualized for sex + SV1-SV3.
#   Age excluded — it IS the signal of interest.
#   β < 0 = repressive; β > 0 = activating.
#
# Input:
#   RNAseq_genes.csv — columns: gene (symbol), cluster (1-4)
#
# Output: /home/tander10/sternerlab/brain_methylation/genebody_GE_stratified/
# ============================================================

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
# Residualization uses ALL 86 samples for stable covariate estimates;
# residuals are subset by age group at the modelling step only.
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
#
# Key difference from promoter script:
#   Promoter script: promoters(upstream=2000, downstream=2000)
#     → TSS−2kb to TSS+2kb only
#   This script: full gene body + 2kb upstream of TSS
#     → TSS−2kb to TES (transcription end site)
#
# Implementation (strand-aware):
#   flank(start=TRUE) gives the 2kb region directly upstream of the TSS,
#   correctly handling strand (upstream of TSS on + strand is to the left;
#   upstream of TSS on − strand is to the right in genomic coordinates).
#   punion() merges the upstream flank with the gene body to give one
#   contiguous region per gene.
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
# flank(start=TRUE) returns the UPSTREAM_KB region immediately 5' of each gene's TSS,
# correctly oriented for both strands. punion() merges it with the gene body.
upstream_flank              <- trim(flank(input_genes_gr, width = UPSTREAM_KB, start = TRUE))
genebody_ext_gr             <- punion(input_genes_gr, upstream_flank, fill.gap = TRUE)
genebody_ext_gr$gene_symbol <- input_genes_gr$symbol

cat(sprintf("  Gene body regions defined: %d genes\n", length(genebody_ext_gr)))
cat(sprintf("  Median region width: %d bp (vs median gene body: %d bp)\n",
            median(width(genebody_ext_gr)),
            median(width(input_genes_gr))))

# ============================================================
# 5. COMPUTE MEAN GENE BODY METHYLATION PER GENE (ALL SAMPLES)
# CpGs within each gene body region are averaged per sample.
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

cat("\nResults saved to:", out_dir, "\n")

# ============================================================
# 12. VISUALIZATIONS
# ============================================================

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

dir_colors  <- c("repressive" = "#2166ac", "activating" = "#d73027")
early_label <- sprintf("Early (age < %d yr)",  AGE_SPLIT)
late_label  <- sprintf("Late  (age >= %d yr)", AGE_SPLIT)

# ---- 12a. Volcano plots ----
make_volcano <- function(results_df, title_str) {
  rp <- results_df[!is.na(results_df$beta_meth), ]
  rp$sig_label <- factor(
    ifelse(rp$sig_05, "FDR < 0.05", ifelse(rp$sig_10, "FDR < 0.10", "NS")),
    levels = c("FDR < 0.05", "FDR < 0.10", "NS")
  )
  ggplot(rp, aes(x = beta_meth, y = -log10(pval_meth),
                  color = direction, alpha = sig_label, size = sig_label)) +
    geom_point() +
    scale_color_manual(values = dir_colors, name = "Direction") +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR < 0.10" = 0.7, "NS" = 0.25),
                       name = "Significance") +
    scale_size_manual(values  = c("FDR < 0.05" = 2, "FDR < 0.10" = 1.5, "NS" = 1),
                      name = "Significance") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      x     = expression(beta ~ "(gene body meth → expression)"),
      y     = expression(-log[10](p-value)),
      title = title_str
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10))
}

p_vol_linear <- make_volcano(results_linear,
  sprintf("Linear genes — all samples (n = %d)", nrow(results_linear)))
p_vol_early  <- make_volcano(results_early,
  sprintf("Nonlinear — %s\n(n = %d genes, %d samples)",
          early_label, nrow(results_early), sum(early_mask)))
p_vol_late   <- make_volcano(results_late,
  sprintf("Nonlinear — %s\n(n = %d genes, %d samples)",
          late_label, nrow(results_late), sum(late_mask)))

ggsave(file.path(fig_dir, "Volcano_Linear_AllSamples.pdf"),
       p_vol_linear, width = 7, height = 5)

p_vol_stratified <- p_vol_early + p_vol_late +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Nonlinear Trajectory Genes: Gene Body Methylation → Expression",
    subtitle = sprintf("Stratified at age = %d years (inflection point)", AGE_SPLIT)
  )
ggsave(file.path(fig_dir, "Volcano_Nonlinear_EarlyVsLate.pdf"),
       p_vol_stratified, width = 12, height = 5)

p_vol_all <- p_vol_linear + p_vol_early + p_vol_late +
  plot_layout(guides = "collect", ncol = 3) +
  plot_annotation(title = "Gene Body Methylation → Expression: All Analyses")
ggsave(file.path(fig_dir, "Volcano_AllAnalyses.pdf"),
       p_vol_all, width = 16, height = 5)

# ---- 12b. Beta distribution histograms ----
make_hist <- function(results_df, title_str) {
  rp <- results_df[!is.na(results_df$beta_meth), ]
  ggplot(rp, aes(x = beta_meth, fill = direction)) +
    geom_histogram(bins = 40, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = dir_colors, name = "Direction") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x        = expression(beta ~ "(gene body meth → expression)"),
      y        = "Number of genes",
      title    = title_str,
      subtitle = sprintf("Repressive: %.1f%%  |  Activating: %.1f%%",
                         100 * mean(rp$direction == "repressive"),
                         100 * mean(rp$direction == "activating"))
    ) +
    theme_bw(base_size = 11)
}

ggsave(file.path(fig_dir, "BetaDist_Linear.pdf"),
       make_hist(results_linear, "Linear genes — all samples"), width = 6, height = 4)
ggsave(file.path(fig_dir, "BetaDist_Nonlinear_Early.pdf"),
       make_hist(results_early, sprintf("Nonlinear — %s", early_label)), width = 6, height = 4)
ggsave(file.path(fig_dir, "BetaDist_Nonlinear_Late.pdf"),
       make_hist(results_late,  sprintf("Nonlinear — %s", late_label)),  width = 6, height = 4)

# ---- 12c. Early vs late beta scatter (nonlinear genes) ----
period_colors <- c(
  "Both periods" = "#7B2D8B",
  "Early only"   = "#1B7837",
  "Late only"    = "#D6604D",
  "Neither"      = "#BDBDBD"
)

comp_df <- comparison %>% dplyr::filter(!is.na(beta_early), !is.na(beta_late))

p_beta_scatter <- ggplot(comp_df,
                          aes(x = beta_early, y = beta_late, color = period_sig)) +
  geom_point(aes(size = period_sig, alpha = period_sig)) +
  scale_color_manual(values = period_colors, name = "Significant at FDR < 0.05") +
  scale_size_manual(
    values = c("Both periods" = 2.5, "Early only" = 2, "Late only" = 2, "Neither" = 1),
    name   = "Significant at FDR < 0.05"
  ) +
  scale_alpha_manual(
    values = c("Both periods" = 1, "Early only" = 0.9, "Late only" = 0.9, "Neither" = 0.3),
    name   = "Significant at FDR < 0.05"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "gray30", linewidth = 0.5) +
  labs(
    x        = expression(beta ~ "early life (age < 10 yr)"),
    y        = expression(beta ~ "late life (age >= 10 yr)"),
    title    = "Gene Body Methylation Effect: Early vs Late Life",
    subtitle = "Nonlinear trajectory genes (Clusters 3+4)\nDotted line = identity (no change across life)"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(fig_dir, "EarlyVsLate_Beta_Scatter.pdf"),
       p_beta_scatter, width = 7, height = 6)

# ---- 12d. Direction summary bars — all three analyses ----
make_dir_bar <- function(results_df, label) {
  rp <- results_df[!is.na(results_df$beta_meth), ]
  tier_rows <- list(
    data.frame(tier = "All tested", direction = "repressive",
               proportion = mean(rp$direction == "repressive")),
    data.frame(tier = "All tested", direction = "activating",
               proportion = mean(rp$direction == "activating")),
    data.frame(tier = "FDR < 0.10", direction = "repressive",
               proportion = if (any(rp$sig_10)) mean(rp$direction[rp$sig_10] == "repressive") else NA),
    data.frame(tier = "FDR < 0.10", direction = "activating",
               proportion = if (any(rp$sig_10)) mean(rp$direction[rp$sig_10] == "activating") else NA),
    data.frame(tier = "FDR < 0.05", direction = "repressive",
               proportion = if (any(rp$sig_05)) mean(rp$direction[rp$sig_05] == "repressive") else NA),
    data.frame(tier = "FDR < 0.05", direction = "activating",
               proportion = if (any(rp$sig_05)) mean(rp$direction[rp$sig_05] == "activating") else NA)
  )
  tier_df <- do.call(rbind, tier_rows) %>%
    dplyr::filter(!is.na(proportion)) %>%
    mutate(tier = factor(tier, levels = c("All tested", "FDR < 0.10", "FDR < 0.05")))

  ggplot(tier_df, aes(x = tier, y = proportion * 100, fill = direction)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = dir_colors, name = "Direction") +
    geom_hline(yintercept = 50, linetype = "dashed",
               color = "white", linewidth = 0.8) +
    labs(x = NULL, y = "% of genes", title = label) +
    theme_bw(base_size = 10)
}

p_dir_all <- make_dir_bar(results_linear, "Linear — all samples") +
             make_dir_bar(results_early,  sprintf("Nonlinear — %s", early_label)) +
             make_dir_bar(results_late,   sprintf("Nonlinear — %s", late_label)) +
  plot_layout(guides = "collect", ncol = 3) +
  plot_annotation(
    title    = "Directional Enrichment of Gene Body Methylation Effects",
    subtitle = "Blue = repressive (higher methylation → lower expression)"
  )
ggsave(file.path(fig_dir, "DirectionSummary_AllAnalyses.pdf"),
       p_dir_all, width = 12, height = 4)

# ---- 12e. Per-gene stratified scatter plots (nonlinear genes) ----
make_stratified_scatter <- function(gene_name, comp_df,
                                    gb_mat, expr_mat, meta,
                                    e_samp, l_samp, age_split) {
  if (!gene_name %in% rownames(gb_mat))   return(NULL)
  if (!gene_name %in% rownames(expr_mat)) return(NULL)

  all_samp <- c(e_samp, l_samp)
  meta_sub <- meta[match(all_samp, meta$sampid), ]

  el <- sprintf("Early (age < %d yr)",  age_split)
  ll <- sprintf("Late  (age >= %d yr)", age_split)

  df <- data.frame(
    meth   = as.numeric(gb_mat[gene_name,   all_samp]),
    expr   = as.numeric(expr_mat[gene_name, all_samp]),
    period = factor(c(rep(el, length(e_samp)), rep(ll, length(l_samp))),
                    levels = c(el, ll)),
    age    = meta_sub$age
  )

  row_info <- comp_df[comp_df$gene == gene_name, ]
  if (nrow(row_info) == 0) return(NULL)

  subtitle_str <- sprintf(
    "Early: β = %.3f, FDR = %.3g  |  Late: β = %.3f, FDR = %.3g",
    row_info$beta_early, row_info$fdr_early,
    row_info$beta_late,  row_info$fdr_late
  )

  period_pal <- setNames(c("#1B7837", "#D6604D"), c(el, ll))

  ggplot(df, aes(x = meth, y = expr, color = period, fill = period, shape = period)) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    scale_color_manual(values = period_pal, name = NULL) +
    scale_fill_manual(values  = period_pal, name = NULL) +
    scale_shape_manual(values = c(16, 17),  name = NULL) +
    labs(
      x        = "Gene body methylation (residual)",
      y        = "Expression (residual)",
      title    = gene_name,
      subtitle = subtitle_str
    ) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle   = element_text(size = 7.5),
          legend.position = "bottom")
}

plot_nonlinear_genes <- unique(c(
  head(early_only$gene, 3),
  head(late_only$gene,  3),
  head(both_sig$gene,   2)
))

if (length(plot_nonlinear_genes) > 0) {
  scatter_list <- Filter(
    Negate(is.null),
    lapply(plot_nonlinear_genes, make_stratified_scatter,
           comp_df   = comparison,
           gb_mat    = gb_nonlinear,
           expr_mat  = expr_nonlinear,
           meta      = metadata_shared,
           e_samp    = early_samples,
           l_samp    = late_samples,
           age_split = AGE_SPLIT)
  )
  if (length(scatter_list) >= 1) {
    ncols <- min(3, length(scatter_list))
    nrows <- ceiling(length(scatter_list) / ncols)
    ggsave(file.path(fig_dir, "Nonlinear_Stratified_Scatter_Examples.pdf"),
           wrap_plots(scatter_list, ncol = ncols),
           width = ncols * 4.5, height = nrows * 4.5)
  }
} else {
  cat("No significant nonlinear genes to plot scatter examples for.\n")
}

# ---- 12f. Top linear gene scatter plots ----
make_linear_scatter <- function(gene_name, results_df, gb_mat, expr_mat, meta) {
  if (!gene_name %in% rownames(gb_mat))   return(NULL)
  if (!gene_name %in% rownames(expr_mat)) return(NULL)

  df <- data.frame(
    meth = as.numeric(gb_mat[gene_name, ]),
    expr = as.numeric(expr_mat[gene_name, ]),
    age  = meta$age
  )
  row_info     <- results_df[results_df$gene == gene_name, ]
  subtitle_str <- sprintf("β = %.3f | FDR = %.3g | %s",
                          row_info$beta_meth, row_info$fdr, row_info$direction)

  ggplot(df, aes(x = meth, y = expr, color = age)) +
    geom_point(size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_viridis_c(option = "plasma", name = "Age (yr)") +
    labs(
      x        = "Gene body methylation (residual)",
      y        = "Expression (residual)",
      title    = gene_name,
      subtitle = subtitle_str
    ) +
    theme_bw(base_size = 11)
}

top_lin_rep  <- head(results_linear$gene[results_linear$direction == "repressive" &
                                           !is.na(results_linear$fdr)][
                       order(results_linear$fdr[results_linear$direction == "repressive" &
                                                  !is.na(results_linear$fdr)])], 3)
top_lin_act  <- head(results_linear$gene[results_linear$direction == "activating" &
                                           !is.na(results_linear$fdr)][
                       order(results_linear$fdr[results_linear$direction == "activating" &
                                                  !is.na(results_linear$fdr)])], 3)
linear_ex    <- unique(c(top_lin_rep, top_lin_act))

if (length(linear_ex) >= 1) {
  lin_list <- Filter(
    Negate(is.null),
    lapply(linear_ex, make_linear_scatter,
           results_df = results_linear,
           gb_mat     = gb_linear,
           expr_mat   = expr_linear,
           meta       = metadata_shared)
  )
  if (length(lin_list) >= 1) {
    ncols <- min(3, length(lin_list))
    ggsave(file.path(fig_dir, "Linear_Scatter_Examples.pdf"),
           wrap_plots(lin_list, ncol = ncols),
           width  = ncols * 4.5,
           height = ceiling(length(lin_list) / ncols) * 4)
  }
}

cat("\nAll figures saved to:", fig_dir, "\n")
cat("Done.\n")
