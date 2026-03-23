# ===============================
# SEX-SPECIFIC DIFFERENTIAL METHYLATION ANALYSIS
# ===============================

library(bsseq)
library(DSS)
library(sva)
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(org.Mmu.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(HDF5Array)

# ---------------------------
# 0. Load preprocessed data
# ---------------------------
output_dir <- "/home/tander10/sternerlab/brain_methylation/results"
bs_file <- "/home/tander10/sternerlab/brain_methylation/results/bs_filtered_86samples_hdf5"
metadata_file <- "/home/tander10/sternerlab/brain_methylation/results/metadata_filtered_86samples_with_SV1-3.rds"

cat("Loading BSseq object (HDF5-backed)...\n")
bs <- loadHDF5SummarizedExperiment(bs_file)

cat("Loading metadata with SVs...\n")
metadata <- readRDS(metadata_file)

# Align metadata to BSseq
stopifnot(all(sampleNames(bs) %in% metadata$sampid))
metadata <- metadata[match(sampleNames(bs), metadata$sampid), ]
stopifnot(all(sampleNames(bs) == metadata$sampid))

# ---------------------------
# 1. Split dataset by Sex
# ---------------------------
sex_groups <- split(metadata$sampid, metadata$Sex)
print(lapply(sex_groups, length))

for (sex in names(sex_groups)) {
  
  cat("\n==============================\n")
  cat(" Running sex-specific model for:", sex, "\n")
  cat("==============================\n\n")
  
  sex_dir <- file.path(output_dir, paste0("SexSpecific_", sex))
  dir.create(sex_dir, showWarnings = FALSE)
  
  # Subset BSseq and metadata
  bs_sex <- bs[, metadata$sampid %in% sex_groups[[sex]]]
  head(rownames(bs_sex))
  meta_sex <- metadata[metadata$sampid %in% sex_groups[[sex]], ]
  rownames(meta_sex) <- meta_sex$sampid
  
  # ---------------------------
  # 2. Run SVA (technical covariates)
  # ---------------------------
  beta_mat <- as.matrix(getMeth(bs_sex, type = "raw"))
  beta_mat[beta_mat == 0] <- 1e-6
  beta_mat[beta_mat == 1] <- 1 - 1e-6
  meth_mat <- log2(beta_mat / (1 - beta_mat))
  meth_mat <- meth_mat[complete.cases(meth_mat), ]
  meth_mat <- meth_mat[apply(meth_mat, 1, var) > 0, ]
  
  cat("Running SVA for", sex, "...\n")
  mod  <- model.matrix(~ Age, data = meta_sex)
  mod0 <- model.matrix(~ 1, data = meta_sex)
  svobj <- sva(meth_mat, mod, mod0)
  
  sv_mat <- svobj$sv
  if (is.null(colnames(sv_mat))) colnames(sv_mat) <- paste0("SV", seq_len(ncol(sv_mat)))
  meta_sex <- cbind(meta_sex, sv_mat)
  
  # Save updated metadata
  saveRDS(meta_sex, file = file.path(sex_dir, paste0("metadata_with_SVs_", sex, ".rds")))
  
  # ---------------------------
  # 3. Build design matrix like variance script
  # ---------------------------
  # Keep up to 3 surrogate variables (if fewer exist, handle gracefully)
  sv_keep <- colnames(sv_mat)[1:min(3, ncol(sv_mat))]
  
  # Build design like your variance analysis
  design_formula <- as.formula(paste("~ Age +", paste(sv_keep, collapse = " + ")))
  design_df <- meta_sex[, c("Age", sv_keep), drop = FALSE]
  
  # Double-check data structure
  stopifnot(is.data.frame(design_df))
  cat("✅ Design matrix columns for", sex, ":", paste(colnames(design_df), collapse = ", "), "\n")
  
  # ---------------------------
  # 4. Run DSS with formula + metadata (not matrix)
  # ---------------------------
  cat("Running DSS for", sex, "...\n")
  DMLfit <- DMLfit.multiFactor(bs_sex, design = design_df, formula = design_formula)
  DMLtest <- DMLtest.multiFactor(DMLfit, coef = "Age")
  saveRDS(DMLtest, file = file.path(sex_dir, paste0("DMLtest_", sex, ".rds")))
  
  # Significant CpGs
  DMLsig <- subset(DMLtest, fdrs < 0.1)
  saveRDS(DMLsig, file = file.path(sex_dir, paste0("DMLsig_", sex, ".rds")))
  cat("Significant CpGs for", sex, ":", nrow(DMLsig), "\n")
  
  # ---------------------------
  # 5. DMR calling
  # ---------------------------
  DMRtest <- callDMR(DMLtest, p.threshold = 0.05, minCG = 3, dis.merge = 300)
  saveRDS(DMRtest, file = file.path(sex_dir, paste0("DMRtest_", sex, ".rds")))
  
  # ---------------------------
  # 6. Annotation
  # ---------------------------
  txdb <- makeTxDbFromEnsembl(organism = "Macaca mulatta", server = "ensembldb.ensembl.org")
  sig_gr <- GRanges(seqnames = DMLsig$chr, ranges = IRanges(DMLsig$pos, DMLsig$pos))
  seqlevelsStyle(sig_gr) <- seqlevelsStyle(txdb)
  
  annotated_sig <- annotatePeak(sig_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
  saveRDS(annotated_sig, file = file.path(sex_dir, paste0("Annotated_DMPs_", sex, ".rds")))
  
  # ---------------------------
  # 7. Manhattan plot
  # ---------------------------
  DMLsig <- DMLsig %>%
    mutate(methylation_status = ifelse(stat > 0, "Hypermethylated", "Hypomethylated"))
  
  manhattan_plot <- ggplot(DMLsig, aes(x = pos, y = log10(abs(stat)), color = methylation_status)) +
    geom_point(alpha = 0.6, size = 0.6) +
    facet_wrap(~chr, scales = "free_x", nrow = 3) +
    scale_color_manual(values = c("Hypermethylated" = "#e34a33", "Hypomethylated" = "#46c19a")) +
    theme_minimal() +
    labs(title = paste("Sex-Specific CpG Sites:", sex),
         x = "Genomic Position", y = "log10 |Test Statistic|") +
    theme(legend.position = "top")
  
  ggsave(file.path(sex_dir, paste0("Manhattan_", sex, ".pdf")), plot = manhattan_plot, width = 10, height = 6)
  
  cat("✅ Completed analysis for", sex, "\n")
  
  # ---------------------------
  # 8. Save significant CpGs as BED files
  # ---------------------------
  if (nrow(DMLsig) > 0) {
    cat("Saving BED files for significant CpGs (", sex, ")...\n")
    
    # Ensure "chr" prefix
    DMLsig$chr <- ifelse(grepl("^chr", DMLsig$chr), DMLsig$chr, paste0("chr", DMLsig$chr))
    
    # Split into hyper and hypo
    hyper_sig <- subset(DMLsig, stat > 0)
    hypo_sig  <- subset(DMLsig, stat < 0)
    
    # --- All significant CpGs ---
    all_bed <- data.frame(
      chrom = DMLsig$chr,
      start = DMLsig$pos,
      end   = DMLsig$pos,
      name  = paste0("CpG_", seq_len(nrow(DMLsig))),
      score = -log10(DMLsig$fdrs),
      strand = "+"
    )
    write.table(all_bed,
                file = file.path(sex_dir, paste0("AllSigCpGs_", sex, ".bed")),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # --- Hypermethylated CpGs ---
    if (nrow(hyper_sig) > 0) {
      hyper_bed <- data.frame(
        chrom = hyper_sig$chr,
        start = hyper_sig$pos,
        end   = hyper_sig$pos,
        name  = paste0("CpG_", seq_len(nrow(hyper_sig))),
        score = -log10(hyper_sig$fdrs),
        strand = "+"
      )
      write.table(hyper_bed,
                  file = file.path(sex_dir, paste0("HypermethylatedCpGs_", sex, ".bed")),
                  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
    # --- Hypomethylated CpGs ---
    if (nrow(hypo_sig) > 0) {
      hypo_bed <- data.frame(
        chrom = hypo_sig$chr,
        start = hypo_sig$pos,
        end   = hypo_sig$pos,
        name  = paste0("CpG_", seq_len(nrow(hypo_sig))),
        score = -log10(hypo_sig$fdrs),
        strand = "+"
      )
      write.table(hypo_bed,
                  file = file.path(sex_dir, paste0("HypomethylatedCpGs_", sex, ".bed")),
                  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
    cat("✅ BED files saved for", sex, "\n")
  } else {
    cat("⚠️ No significant CpGs to export for", sex, "\n")
  }
}

cat("\n✅ All sex-specific models completed successfully!\n")
