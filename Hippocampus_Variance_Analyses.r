#!/usr/bin/env Rscript
# ==============================
# Variance Analysis with missMethyl (with SV1–SV3 adjustment)
# ==============================

# ---- Libraries ----
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
library(txdbmaker)

# ---- Paths ----
output_dir <- "/home/tander10/sternerlab/brain_methylation/results_variance"
dir.create(output_dir, showWarnings = FALSE)

bs_file <- "/home/tander10/sternerlab/brain_methylation/results/bs_filtered_86samples_hdf5"
metadata_file <- "/home/tander10/sternerlab/brain_methylation/results/metadata_filtered_86samples_with_SV1-3.rds"

# ==== Step 1: Load Data ====
cat("Loading BSseq object (HDF5-backed)...\n")
bs <- loadHDF5SummarizedExperiment(bs_file)

cat("Loading metadata with SVs...\n")
metadata <- readRDS(metadata_file)

# ==== Step 2: Ensure alignment ====
if (!"sampid" %in% colnames(metadata)) stop("Metadata must contain 'sampid' column.")
metadata <- metadata[match(sampleNames(bs), metadata$sampid), ]
stopifnot(all(sampleNames(bs) == metadata$sampid))

required_cols <- c("Age", "Sex", "SV1", "SV2", "SV3")
stopifnot(all(required_cols %in% colnames(metadata)))

# ==== Step 3: Prepare M-values ====
cat("Extracting M-values...\n")
meth <- as.matrix(bsseq::getMeth(bs, type = "raw"))
rownames(meth) <- paste0(as.character(seqnames(bs)), "_", start(bs), "_", end(bs))
meth[meth == 0] <- 1e-6
meth[meth == 1] <- 1 - 1e-6
meth_M <- log2(meth / (1 - meth))
meth_M <- meth_M[complete.cases(meth_M), ]

# ==== Step 4: Design matrix ====
metadata$Sex <- as.factor(metadata$Sex)
design_var <- model.matrix(~ Age + Sex + SV1 + SV2 + SV3, data = metadata)

# ==== Step 5: Run variance analysis ====
cat("Running variance analysis...\n")
varfit <- missMethyl::varFit(meth_M, design = design_var, coef = "Age", trend = TRUE)
allvar_results <- missMethyl::topVar(varfit, coef = "Age", number = nrow(meth_M))
allvar_results$CpG_ID <- rownames(allvar_results)
allvar_results$Age_FDR <- p.adjust(allvar_results$P.Value, method = "fdr")

filtered_results <- subset(allvar_results, Age_FDR < 0.1)
cat("Significant CpGs (FDR < 0.1):", nrow(filtered_results), "\n")

# ==== Step 6: Save core results ====
saveRDS(allvar_results, file.path(output_dir, "topVar_Age_Sex_SV1-3_all_results.rds"))
write.csv(allvar_results, file.path(output_dir, "topVar_Age_Sex_SV1-3_all_results.csv"), row.names = FALSE)
saveRDS(filtered_results, file.path(output_dir, "topVar_Age_Sex_SV1-3_FDR0.1_results.rds"))
write.csv(filtered_results, file.path(output_dir, "topVar_Age_Sex_SV1-3_FDR0.1_results.csv"), row.names = FALSE)

# ==== Step 7: Stratify and Save Increasing/Decreasing + All ====
increasing_var <- filtered_results[filtered_results$LogVarRatio > 0, ]
decreasing_var <- filtered_results[filtered_results$LogVarRatio < 0, ]

cat("VMPs with increasing variance:", nrow(increasing_var), "\n")
cat("VMPs with decreasing variance:", nrow(decreasing_var), "\n")

# Save all, increasing, and decreasing VMPs as CSV
write.csv(filtered_results, file.path(output_dir, "AllVMPs_Age.csv"), row.names = FALSE)
write.csv(increasing_var, file.path(output_dir, "VMPs_IncreasingVariance_Age.csv"), row.names = FALSE)
write.csv(decreasing_var, file.path(output_dir, "VMPs_DecreasingVariance_Age.csv"), row.names = FALSE)

# ==== Step 8: Export BED files ====
cpg_to_bed <- function(cpg_ids) {
  bed <- do.call(rbind, strsplit(cpg_ids, "_"))
  bed <- as.data.frame(bed, stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end")
  bed$chr <- paste0("chr", bed$chr)
  bed
}

if (nrow(filtered_results) > 0) {
  bed_all <- cpg_to_bed(filtered_results$CpG_ID)
  write.table(bed_all, file.path(output_dir, "AllVMPs_Age.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

if (nrow(increasing_var) > 0) {
  bed_inc <- cpg_to_bed(increasing_var$CpG_ID)
  write.table(bed_inc, file.path(output_dir, "VMPs_IncreasingVariance_Age.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

if (nrow(decreasing_var) > 0) {
  bed_dec <- cpg_to_bed(decreasing_var$CpG_ID)
  write.table(bed_dec, file.path(output_dir, "VMPs_DecreasingVariance_Age.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

cat("✅ BED files written for all, increasing, and decreasing VMPs.\n")

# ==== Step 9: Annotation of all VMPs ====
if (nrow(filtered_results) > 0) {
  cat("🔍 Annotating all VMPs using ChIPseeker...\n")
  
  txdb <- txdbmaker::makeTxDbFromEnsembl(
    organism = "Macaca mulatta",
    server = "ensembldb.ensembl.org"
  )
  
  # Convert all VMPs to GRanges
  bed_all <- cpg_to_bed(filtered_results$CpG_ID)
  gr_all <- GRanges(seqnames = bed_all$chr,
                    ranges = IRanges(start = as.numeric(bed_all$start),
                                     end = as.numeric(bed_all$end)))
  
  # Harmonize seqlevels safely
  seqlevelsStyle(gr_all) <- seqlevelsStyle(txdb)
  common_levels <- intersect(seqlevels(gr_all), seqlevels(txdb))
  gr_all <- keepSeqlevels(gr_all, common_levels, pruning.mode = "coarse")
  
  # Annotate all
  annot_all <- annotatePeak(
    gr_all,
    TxDb = txdb,
    tssRegion = c(-3000, 3000),
    annoDb = "org.Mmu.eg.db"
  )
  
  saveRDS(annot_all, file.path(output_dir, "Annotated_AllVMPs_Age.rds"))
  cat("✅ Annotation completed for all VMPs.\n")
}

# ==== Step 10: Annotation for Increasing/Decreasing (same style) ====
if (nrow(increasing_var) > 0 | nrow(decreasing_var) > 0) {
  cat("🔍 Annotating increasing/decreasing variance CpGs...\n")
  
  cpg_to_gr <- function(df) {
    split <- do.call(rbind, strsplit(df$CpG_ID, "_"))
    GRanges(seqnames = paste0("chr", split[, 1]),
            ranges = IRanges(start = as.numeric(split[, 2]),
                             end = as.numeric(split[, 3])))
  }
  
  if (nrow(increasing_var) > 0) {
    inc_gr <- cpg_to_gr(increasing_var)
    seqlevelsStyle(inc_gr) <- seqlevelsStyle(txdb)
    inc_gr <- keepSeqlevels(inc_gr, intersect(seqlevels(inc_gr), seqlevels(txdb)), pruning.mode = "coarse")
    annot_inc <- annotatePeak(inc_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
    saveRDS(annot_inc, file.path(output_dir, "Annotated_VMPs_IncreasingVariance.rds"))
  }
  
  if (nrow(decreasing_var) > 0) {
    dec_gr <- cpg_to_gr(decreasing_var)
    seqlevelsStyle(dec_gr) <- seqlevelsStyle(txdb)
    dec_gr <- keepSeqlevels(dec_gr, intersect(seqlevels(dec_gr), seqlevels(txdb)), pruning.mode = "coarse")
    annot_dec <- annotatePeak(dec_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
    saveRDS(annot_dec, file.path(output_dir, "Annotated_VMPs_DecreasingVariance.rds"))
  }
  
  cat("✅ ChIPseeker annotation completed for variance CpGs.\n")
}

cat("✅ All files (All/Increasing/Decreasing) saved in:", output_dir, "\n")
