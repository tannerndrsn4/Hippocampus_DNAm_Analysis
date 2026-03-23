#!/usr/bin/env Rscript
# ===============================
# SEX-SPECIFIC VARIANCE ANALYSIS (missMethyl)
# ===============================

library(bsseq)
library(missMethyl)
library(comethyl)
library(HDF5Array)
library(ChIPseeker)
library(org.Mmu.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)

# ---------------------------
# 0. Paths and setup
# ---------------------------
output_dir <- "/home/tander10/sternerlab/brain_methylation/results_variance_sex"
dir.create(output_dir, showWarnings = FALSE)

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
  cat(" Running sex-specific variance model for:", sex, "\n")
  cat("==============================\n\n")
  
  sex_dir <- file.path(output_dir, paste0("Variance_", sex))
  dir.create(sex_dir, showWarnings = FALSE)
  
  # Subset BSseq and metadata
  bs_sex <- bs[, metadata$sampid %in% sex_groups[[sex]]]
  meta_sex <- metadata[metadata$sampid %in% sex_groups[[sex]], ]
  rownames(meta_sex) <- meta_sex$sampid
  
  # ---------------------------
  # 2. Prepare M-values
  # ---------------------------
  cat("Extracting M-values for", sex, "...\n")
  meth <- as.matrix(bsseq::getMeth(bs_sex, type = "raw"))
  rownames(meth) <- paste0(as.character(seqnames(bs_sex)), "_", start(bs_sex), "_", end(bs_sex))
  meth[meth == 0] <- 1e-6
  meth[meth == 1] <- 1 - 1e-6
  meth_M <- log2(meth / (1 - meth))
  meth_M <- meth_M[complete.cases(meth_M), ]
  
  # ---------------------------
  # 3. Build design matrix (Age + up to 3 SVs)
  # ---------------------------
  sv_keep <- intersect(c("SV1", "SV2", "SV3"), colnames(meta_sex))
  design_var <- model.matrix(as.formula(paste("~ Age +", paste(sv_keep, collapse = " + "))), data = meta_sex)
  
  # ---------------------------
  # 4. Run variance analysis
  # ---------------------------
  cat("Running variance analysis for", sex, "...\n")
  varfit <- varFit(meth_M, design = design_var, coef = "Age", trend = TRUE)
  allvar_results <- topVar(varfit, coef = "Age", number = nrow(meth_M))
  
  # Add IDs and FDR
  allvar_results$CpG_ID <- rownames(allvar_results)
  allvar_results$Age_FDR <- p.adjust(allvar_results$P.Value, method = "fdr")
  
  filtered_results <- subset(allvar_results, Age_FDR < 0.1)
  cat("Significant CpGs (FDR < 0.1) for", sex, ":", nrow(filtered_results), "\n")
  
  # ---------------------------
  # 5. Stratify by variance direction
  # ---------------------------
  increasing_var <- subset(filtered_results, LogVarRatio > 0)
  decreasing_var <- subset(filtered_results, LogVarRatio < 0)
  
  # ---------------------------
  # 6. Save results (all / inc / dec)
  # ---------------------------
  saveRDS(filtered_results, file = file.path(sex_dir, paste0("AllVMPs_", sex, ".rds")))
  write.csv(filtered_results, file = file.path(sex_dir, paste0("AllVMPs_", sex, ".csv")), row.names = FALSE)
  
  write.csv(increasing_var, file.path(sex_dir, paste0("VMPs_IncreasingVariance_Age_", sex, ".csv")), row.names = FALSE)
  write.csv(decreasing_var, file.path(sex_dir, paste0("VMPs_DecreasingVariance_Age_", sex, ".csv")), row.names = FALSE)
  
  # BED helper
  cpg_to_bed <- function(cpg_ids) {
    bed <- do.call(rbind, strsplit(cpg_ids, "_"))
    bed <- as.data.frame(bed, stringsAsFactors = FALSE)
    colnames(bed) <- c("chr", "start", "end")
    bed$chr <- paste0("chr", bed$chr)
    bed
  }
  
  if (nrow(filtered_results) > 0) {
    all_bed <- cpg_to_bed(filtered_results$CpG_ID)
    write.table(all_bed, file.path(sex_dir, paste0("AllVMPs_", sex, ".bed")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if (nrow(increasing_var) > 0) {
    bed_inc <- cpg_to_bed(increasing_var$CpG_ID)
    write.table(bed_inc, file.path(sex_dir, paste0("VMPs_IncreasingVariance_Age_", sex, ".bed")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if (nrow(decreasing_var) > 0) {
    bed_dec <- cpg_to_bed(decreasing_var$CpG_ID)
    write.table(bed_dec, file.path(sex_dir, paste0("VMPs_DecreasingVariance_Age_", sex, ".bed")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # ---------------------------
  # 7. Annotation (safe version, no seqlevel remapping)
  # ---------------------------
  # ---- Annotation (fixed seqlevel mapping) ----
  if (nrow(filtered_results) > 0) {
    cat("Annotating VMPs for", sex, "...\n")
    
    # Build TxDb for macaque (Ensembl)
    txdb <- txdbmaker::makeTxDbFromEnsembl(
      organism = "Macaca mulatta",
      server = "ensembldb.ensembl.org"
    )
    
    # Convert CpG coordinates into GRanges
    all_bed <- cpg_to_bed(filtered_results$CpG_ID)
    gr <- GRanges(
      seqnames = all_bed$chr,
      ranges = IRanges(start = as.numeric(all_bed$start),
                       end   = as.numeric(all_bed$end))
    )
    
    # Align naming styles between your GRanges and TxDb
    seqlevelsStyle(gr) <- seqlevelsStyle(txdb)
    
    # Drop any seqlevels not found in TxDb to prevent merge errors
    common_levels <- intersect(seqlevels(gr), seqlevels(txdb))
    gr <- keepSeqlevels(gr, common_levels, pruning.mode = "coarse")
    
    # Annotate peaks safely
    annotated_vmp <- annotatePeak(
      gr,
      TxDb = txdb,
      tssRegion = c(-3000, 3000),
      annoDb = "org.Mmu.eg.db"
    )
    
    saveRDS(annotated_vmp, file = file.path(sex_dir, paste0("Annotated_AllVMPs_", sex, ".rds")))
    cat("✅ Annotation completed for", sex, "\n")
  }
  # ---------------------------
  # 8. Visualization
  # ---------------------------
  if (nrow(filtered_results) > 0) {
    p1 <- ggplot(filtered_results, aes(x = LogVarRatio)) +
      geom_histogram(bins = 80, fill = "#1f78b4", alpha = 0.7) +
      theme_minimal() +
      labs(title = paste0("Distribution of LogVarRatio (", sex, ")"),
           x = "LogVarRatio", y = "Count")
    ggsave(file.path(sex_dir, paste0("LogVarRatio_Histogram_", sex, ".pdf")), plot = p1, width = 7, height = 5)
  }
  
  cat("✅ Completed variance analysis for", sex, "\n")
}

cat("\n✅ All sex-specific variance models completed successfully!\n")
