# Load Required Libraries
library(Biostrings)
library(bsseq)
library(comethyl)
library(DSS)
library(ChIPseeker)
library(AnnotationHub)
library(GenomicFeatures)
library(org.Mmu.eg.db)
library(GenomicRanges)
library(Gviz)
library(txdbmaker)
library(RMariaDB)
library(ggplot2)
library(gridExtra)

# Define Paths
output_dir <- "/home/tander10/sternerlab/brain_methylation/results"
dir.create(output_dir, showWarnings = FALSE)

# Step 1: Data Loading and Preprocessing
mmul10_fa <- readRDS(file = "/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")
names(mmul10_fa)[1:21] -> chr
mmul10_fa <- mmul10_fa[chr]
loci <- findLoci("CG", subject = mmul10_fa)

bismarkBSseq <- read.bismark(
  files = list.files("/home/tander10/sternerlab/brain_methylation/Bismark_output",
                     full.names = TRUE, pattern = ".cov.gz"),
  strandCollapse = FALSE, verbose = TRUE, BACKEND = "HDF5Array", loci = loci,
  rmZeroCov = FALSE, dir = "/projects/sternerlab/shared/Coverage_files/hdf5",
  replace = TRUE, BPPARAM = BiocParallel::SerialParam()
)

# Step 2: Filtering
bs <- filterCpGs(
  bismarkBSseq, cov = 10, perSample = 0.8,
  file = file.path(output_dir, "Filtered_BSseq_Hippocampus_5X.rds")
)
#bs <- readRDS("/home/tander10/sternerlab/brain_methylation/results/Filtered_BSseq_Hippocampus_5X.rds")
#sampleNames(bs)

# ----------------------
# Step 2.5: Relabel and Subset by sampid
# ----------------------

# Load metadata (original = 96 samples, with SubjectID instead of sampid)
meta_original <- read.csv("/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata.csv")

# Rename SubjectID column to sampid for consistency
if ("SubjectID" %in% colnames(meta_original) && !"sampid" %in% colnames(meta_original)) {
  colnames(meta_original)[colnames(meta_original) == "SubjectID"] <- "sampid"
}

# 1. Relabel bsseq sample names to sampid values
if (length(sampleNames(bs)) != nrow(meta_original)) {
  stop("Number of samples in bsseq object does not match metadata rows. 
       Check ordering before renaming.")
}

sampleNames(bs) <- meta_original$sampid
cat("✅ Renamed bsseq sampleNames to sampid values.\n")
print(head(sampleNames(bs)))

# 2. Load RNA-seq metadata
rnaseq_meta <- readRDS("/home/tander10/sternerlab/brain_methylation/RNAseq_metadata.rds")
keep_ids <- unique(rnaseq_meta$sampid)

# 3. Subset bsseq and metadata to only sampids present in RNA-seq metadata
common_samples <- intersect(sampleNames(bs), keep_ids)
bs <- bs[, common_samples]

metadata <- meta_original[meta_original$sampid %in% common_samples, ]
metadata <- droplevels(metadata)

# 4. Verify alignment
stopifnot(all(sampleNames(bs) == metadata$sampid))

cat("✅ Number of samples retained after subsetting:", nrow(metadata), "\n")

# Simple stacked histogram of Age by Sex
samp <- ggplot(metadata, aes(x = Age, fill = Sex)) +
  geom_histogram(binwidth = 5, position = "stack", color = "black") +
  labs(
    x = "Age (years)",
    y = "Frequency of adult samples"
  ) +
  scale_fill_manual(values = c("M" = "#6980ce", "F" = "#b561b8")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Save the plot
ggsave(
  filename = file.path(output_dir, "Age_Distribution_by_Sex.pdf"),
  plot = samp,
  width = 8, height = 5, units = "in"
)

# Save aligned bs and metadata for downstream analyses (e.g., variance analysis)
library(HDF5Array)

saveHDF5SummarizedExperiment(
  bs,
  dir = file.path(output_dir, "bs_filtered_86samples_hdf5"),
  replace = TRUE
)

# ==========================
# Coverage QC plot after filtering
# ==========================
cat("Creating CpG coverage distribution plot...\n")

# Extract coverage matrix: rows = CpGs, columns = samples
cov_matrix <- getCoverage(bs, type = "Cov")

# Convert to long format for plotting
cov_df <- as.data.frame(as.matrix(cov_matrix))
cov_df$CpG <- rownames(cov_df)
cov_long <- tidyr::pivot_longer(
  cov_df,
  cols = -CpG,
  names_to = "SampleID",
  values_to = "Coverage"
)

# Optional: log-transform coverage for visualization
cov_long$logCov <- log10(cov_long$Coverage + 1)

# Plot distribution of coverage per sample
coverage_plot <- ggplot(cov_long, aes(x = SampleID, y = logCov)) +
  geom_violin(fill = "steelblue", alpha = 0.6, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "CpG Coverage Distribution per Sample (after filtering)",
    x = "Sample",
    y = "log10(Coverage + 1)"
  )

# Save plot
ggsave(
  filename = file.path(output_dir, "CpG_Coverage_Distribution_After_Filtering.pdf"),
  plot = coverage_plot,
  width = 12, height = 6, units = "in"
)

### Section for identifying and controlling for technical covariates that do not correlate with age or sex ###

library(sva)

# Step 1 — Ensure filtering matches DSS
# (bs object here is already filtered with cov = 5, perSample = 0.6 in Step 2)

# Step 2 — Extract beta values (0–1) for these CpGs
beta_mat <- as.matrix(getMeth(bs, type = "raw"))

# Bound values away from 0/1 to avoid infinities in logit transform
beta_mat[beta_mat == 0] <- 1e-6
beta_mat[beta_mat == 1] <- 1 - 1e-6

# Step 3 — Convert to M-values for statistical modeling
meth_mat <- log2(beta_mat / (1 - beta_mat))

# Step 4 — Filter sites for SVA
meth_mat <- meth_mat[complete.cases(meth_mat), ]                       # remove missing
dim(meth_mat)
meth_mat <- meth_mat[apply(meth_mat, 1, var) > 0, ]                    # remove constant sites

# Keep top variable CpGs (optional for speed; adjust number if needed)
top_n <- 10000
if (nrow(meth_mat) > top_n) {
  top_idx <- order(apply(meth_mat, 1, var), decreasing = TRUE)[1:top_n]
  meth_mat <- meth_mat[top_idx, ]
}

# Step 5 — Run SVA with Age + Sex in the model
mod  <- model.matrix(~ Age + Sex, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)

# Use sva()
svobj <- sva(meth_mat, mod, mod0)

# Step 6 — Extract and name surrogate variables
sv_mat <- svobj$sv
if (is.null(colnames(sv_mat))) {
  colnames(sv_mat) <- paste0("SV", seq_len(ncol(sv_mat)))
}

# Step 7 — Correlation testing with Age and Sex
metadata_num <- metadata
metadata_num$Sex <- as.numeric(as.factor(metadata_num$Sex))
metadata_num$Age <- as.numeric(metadata_num$Age)

safe_cor <- function(x, y) {
  if (all(is.na(x)) || all(is.na(y))) return(NA)
  if (var(x, na.rm = TRUE) == 0 || var(y, na.rm = TRUE) == 0) return(NA)
  suppressWarnings(cor(x, y, use = "complete.obs"))
}

sv_assoc <- data.frame(SV = colnames(sv_mat))
sv_assoc$Age_cor <- apply(sv_mat, 2, function(sv) cor(metadata_num$Age, sv, use = "complete.obs"))
sv_assoc$Age_p   <- apply(sv_mat, 2, function(sv) cor.test(metadata_num$Age, sv)$p.value)
sv_assoc$Sex_p   <- apply(sv_mat, 2, function(sv) t.test(sv ~ metadata$Sex)$p.value)

sv_assoc$Age_FDR <- p.adjust(sv_assoc$Age_p, method = "fdr")
sv_assoc$Sex_FDR <- p.adjust(sv_assoc$Sex_p, method = "fdr")

# Step 8 — Keep only SVs that are not significantly associated with Age/Sex
sv_technical <- sv_assoc$SV[sv_assoc$Age_FDR >= 0.05 & sv_assoc$Sex_FDR >= 0.05]

# Step 9 — Keep only high-variance SVs (adjust threshold if needed)
sv_var  <- apply(sv_mat, 2, var)
sv_prop <- sv_var / sum(sv_var)
threshold <- 0.02  # 2% of SV variance
sv_highvar <- sv_technical[sv_prop[sv_technical] > threshold]

# Calculate variance and proportion for each SV
sv_var  <- apply(sv_mat, 2, var)
sv_prop <- sv_var / sum(sv_var)

# Create variance plot
sv_df <- data.frame(
  SV = seq_along(sv_prop),
  CumulativeVariance = cumsum(sv_prop)
)

# Save the plot to your output directory
pdf(file.path(output_dir, "SV_cumulative_variance.pdf"), width = 6, height = 5)
plot(
  sv_df$SV,
  sv_df$CumulativeVariance,
  type = "b",
  xlab = "Number of SVs",
  ylab = "Cumulative Proportion of Variance Explained",
  main = "Cumulative Variance Explained by SVs"
)
abline(h = 0.8, col = "red", lty = 2)  # Example threshold
dev.off()

# Calculate R² for each SV against the methylation matrix
sv_r2 <- apply(sv_mat, 2, function(sv) {
  # Correlation-based R² between SV and each CpG site, averaged across sites
  # This avoids lm() on a giant matrix, which would be slow
  cor_sq <- apply(meth_mat, 1, function(row) {
    suppressWarnings(cor(row, sv, use = "complete.obs")^2)
  })
  mean(cor_sq, na.rm = TRUE)
})

# Save R² barplot to PDF
pdf(file.path(output_dir, "SV_variance_barplot.pdf"), width = 8, height = 5)
barplot(
  sv_r2,
  names.arg = paste0("SV", seq_along(sv_r2)),
  xlab = "Surrogate Variable",
  ylab = "Proportion of Variance Explained in Methylation Data)",
  main = "Variance Explained by Each SV",
  las = 2,               # Rotate x-axis labels
  col = "steelblue"
)
abline(h = 0.02, col = "red", lty = 2)  # Example 5% threshold
dev.off()

cat("SVs to control for:", sv_highvar, "\n")

# Step 10 — Add selected SVs to metadata for DSS modeling
metadata <- cbind(metadata, sv_mat[, sv_highvar, drop = FALSE])

# Optional: QC — see how SVs relate to top PCs
pc <- prcomp(t(meth_mat), scale. = TRUE)
pc_scores <- pc$x[, 1:10]
sv_pc_cor <- cor(sv_mat, pc_scores)
heatmap(sv_pc_cor, main = "Correlation Between SVs and Top PCs")



# Remove samples (columns) with zero variance
# Extract methylation matrix
meth_mat <- as.matrix(getMeth(bs, type = "raw"))

# Keep CpGs with no missing values
meth_mat <- meth_mat[complete.cases(meth_mat), ]
dim(meth_mat)
saveRDS(meth_mat, file = file.path(output_dir, "meth_matrix.rds"))

# Calculate variance for each CpG across samples
cpg_vars <- apply(meth_mat, 1, var, na.rm = TRUE)

# Keep top 10,000 most variable CpGs
top_n <- 10000
if (nrow(meth_mat) > top_n) {
  top_cpgs <- order(cpg_vars, decreasing = TRUE)[1:top_n]
  meth_mat <- meth_mat[top_cpgs, , drop = FALSE]
}

# Remove samples (columns) with zero variance
non_constant_cols <- apply(meth_mat, 2, function(x) var(x, na.rm = TRUE) > 0)
meth_mat <- meth_mat[, non_constant_cols, drop = FALSE]

# PCA
pc <- prcomp(t(meth_mat), scale. = TRUE)
pc_scores <- pc$x[, 1:10]

sv_pc_cor <- cor(sv_mat, pc_scores)
sv_pc_cor
heatmap(sv_pc_cor, main = "Correlation Between SVs and Top PCs")
###

# Step 3: Differential Methylation Analysis - Site-Based

# --- Select SVs to control for ---
sv_keep <- paste0("SV", 1:3)  # SV1, SV2, SV3

# Add them to metadata (if not already present)
metadata <- cbind(metadata, sv_mat[, sv_keep, drop = FALSE])
saveRDS(metadata, file = file.path(output_dir, "metadata_with_SV1-3.rds"))

# Build design matrix including SV1–SV3
design_matrix <- model.matrix(~ Age + Sex + SV1 + SV2 + SV3, data = metadata)
design_matrix <- as.data.frame(design_matrix)

# Run DSS with selected SVs
DMLfit <- DMLfit.multiFactor(bs, design = design_matrix, formula = ~ Age + SexM + SV1 + SV2 + SV3)
DMLtest <- DMLtest.multiFactor(DMLfit, coef = "Age")

# Save results
saveRDS(DMLtest, file = file.path(output_dir, "DMLtest_results_with_SV1-3.rds"))

# ✅ Keep only significant sites (FDR < 0.1) for Manhattan plot
DMLtest_significant <- subset(DMLtest, fdrs < 0.1)
summary(DMLtest$fdrs < 0.1)

# ✅ Split into hyper- and hypomethylated sites
hyper_dms <- subset(DMLtest_significant, stat > 0)
hypo_dms <- subset(DMLtest_significant, stat < 0)

# Create GRanges objects
hyper_dms_gr <- GRanges(seqnames = hyper_dms$chr, ranges = IRanges(hyper_dms$pos, hyper_dms$pos))
hypo_dms_gr <- GRanges(seqnames = hypo_dms$chr, ranges = IRanges(hypo_dms$pos, hypo_dms$pos))

# Step 4: Differential Methylation Analysis - Region-Based
DMRtest_ARIMA <- callDMR(DMLtest, p.threshold = 0.05, minCG = 3, dis.merge = 300)
saveRDS(DMRtest_ARIMA, file = file.path(output_dir, "DMRtest_ARIMA_results.rds"))
dim(DMRtest_ARIMA)

DMRtest <- callDMR(DMLtest_significant, p.threshold = 0.05, minCG = 3, dis.merge = 300)
saveRDS(DMRtest, file = file.path(output_dir, "DMRtest_results.rds"))
dim(DMRtest)

# ✅ Split into hyper- and hypomethylated DMRs
hyper_dmrs <- subset(DMRtest, areaStat > 0)
hypo_dmrs <- subset(DMRtest, areaStat < 0)

# Create GRanges objects
hyper_dmrs_gr <- GRanges(seqnames = hyper_dmrs$chr, ranges = IRanges(hyper_dmrs$start, hyper_dmrs$end))
hypo_dmrs_gr <- GRanges(seqnames = hypo_dmrs$chr, ranges = IRanges(hypo_dmrs$start, hypo_dmrs$end))

# Step 5: Annotation with ChIPseeker
txdb <- makeTxDbFromEnsembl(organism = "Macaca mulatta", server = "ensembldb.ensembl.org")
#txdb <- loadDb("Macaca_mulatta_txdb.sqlite")

# ✅ Annotate background CpGs (All tested sites)
background_cpgs_gr <- GRanges(seqnames = DMLtest$chr, ranges = IRanges(DMLtest$pos, DMLtest$pos))
saveRDS(background_cpgs_gr, file = file.path(output_dir, "annotated_background.rds"))

# Annotate sites and DMRs
annotated_hyper_dms <- annotatePeak(hyper_dms_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_hypo_dms <- annotatePeak(hypo_dms_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_hyper_dmrs <- annotatePeak(hyper_dmrs_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_hypo_dmrs <- annotatePeak(hypo_dmrs_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_background <- annotatePeak(background_cpgs_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")

saveRDS(annotated_hyper_dms, file = file.path(output_dir, "hyper_annotated_cpgs.rds"))
saveRDS(annotated_hypo_dms, file = file.path(output_dir, "hypo_annotated_cpgs.rds"))
saveRDS(annotated_hyper_dmrs, file = file.path(output_dir, "hyper_annotated_dmrs.rds"))
saveRDS(annotated_hypo_dmrs, file = file.path(output_dir, "hypo_annotated_dmrs.rds"))
saveRDS(annotated_background, file = file.path(output_dir, "annotated_background_DMR.rds"))

# Step 6: Visualization

# ✅ Manhattan Plot (Only Significant Sites, Faceted by Chromosome)
manhattan_sites <- ggplot(DMLtest_significant, aes(x = pos, y = log10(abs(stat)), color = stat > 0)) +
  geom_point(alpha = 0.6, size = 0.5) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +  # Separate panels per chromosome
  scale_color_manual(values = c("FALSE" = "#DC3220", "TRUE" = "#005AB5"),
                     labels = c("Hypomethylated", "Hypermethylated")) +
  theme_minimal() +
  labs(
    title = "Manhattan Plot of Significant CpG Sites",
    x = "Genomic Position",
    y = "log10 |Test Statistic|",
    color = "Methylation Status"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 10)
  )

ggsave(file.path(output_dir, "Manhattan_Plot_Sites_Significant.pdf"), plot = manhattan_sites)

# ✅ Manhattan Plot (Compact Faceted by Chromosome)
manhattan_sites <- ggplot(DMLtest_significant, aes(x = pos, y = log10(abs(stat)), color = stat > 0)) +
  geom_point(alpha = 0.6, size = 0.5) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +  # Removed space = "free_x"
  scale_color_manual(values = c("FALSE" = "#DC3220", "TRUE" = "#005AB5"),
                     labels = c("Hypomethylated", "Hypermethylated")) +
  theme_minimal() +
  labs(
    title = "Manhattan Plot of Significant CpG Sites",
    x = "Genomic Position",
    y = "log10 |Test Statistic|",
    color = "Methylation Status"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 9),
    panel.spacing.x = unit(0.1, "lines")
  )

# ✅ Save with smaller width
ggsave(
  file.path(output_dir, "Manhattan_Plot_Sites_Significant.pdf"),
  plot = manhattan_sites,
  width = 12, height = 6, units = "in"  # Shrink width to fit more panels
)

library(ggplot2)
library(dplyr)

# Prepare data: ensure methylation status exists
DMLtest_significant <- DMLtest_significant %>%
  mutate(
    methylation_status = ifelse(stat > 0, "Hypermethylated", "Hypomethylated"),
    chr = factor(chr, levels = sort(unique(chr)))
  )

# Summarize counts of hypo/hyper per chromosome
bar_data <- DMLtest_significant %>%
  group_by(chr, methylation_status) %>%
  tally(name = "count")

# Create bar plot (raw counts instead of proportions)
bar_plot <- ggplot(bar_data, aes(x = chr, y = count, fill = methylation_status)) +
  geom_bar(stat = "identity", position = "stack") +  # stack shows raw counts
  scale_fill_manual(values = c("Hypermethylated" = "#e34a33", "Hypomethylated" = "#46c19a")) +
  theme_minimal() +
  labs(
    x = "Chromosome",
    y = "Count",
    fill = "Methylation Status"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save updated plot
ggsave(
  file.path(output_dir, "Bar_Plot_Hyper_Hypo_Counts.pdf"),
  bar_plot,
  width = 10, height = 5, units = "in"
)

# Prepare data: ensure methylation status exists
DMLtest_significant <- DMLtest_significant %>%
  mutate(
    methylation_status = ifelse(stat > 0, "Hypermethylated", "Hypomethylated"),
    chr = factor(chr, levels = sort(unique(chr)))
  )

# Summarize counts of hypo/hyper per chromosome
bar_data <- DMLtest_significant %>%
  group_by(chr, methylation_status) %>%
  tally(name = "count")

# Create bar plot (single panel)
bar_plot <- ggplot(bar_data, aes(x = chr, y = count, fill = methylation_status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Hypermethylated" = "#b54d88", "Hypomethylated" = "#46c19a")) +
  theme_minimal() +
  labs(
    title = "Proportion of Hyper- vs Hypomethylated CpG Sites per Chromosome",
    x = "Chromosome",
    y = "Proportion",
    fill = "Methylation Status"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plot
ggsave(
  file.path(output_dir, "Bar_Plot_Hyper_Hypo_Proportion.pdf"),
  bar_plot,
  width = 10, height = 5, units = "in"
)


# ✅ Manhattan Plot for DMRs
manhattan_dmrs <- ggplot(DMRtest, aes(x = start, y = abs(areaStat), color = factor(chr))) +
  geom_point(alpha = 0.6) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +
  theme_minimal() +
  labs(title = "Manhattan Plot of DMRs", x = "Genomic Position", y = "Effect Size (|AreaStat|)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(output_dir, "Manhattan_Plot_DMRs.pdf"), plot = manhattan_dmrs)

# ✅ Annotation Bar Plot (Including Background CpGs)
anno_combined <- plotAnnoBar(
  list(
    "Hyper DMPs" = annotated_hyper_dms,
    "Hypo DMPs" = annotated_hypo_dms,
    "Hyper DMRs" = annotated_hyper_dmrs,
    "Hypo DMRs" = annotated_hypo_dmrs,
    "Background CpGs" = annotated_background
  )
)
ggsave(file.path(output_dir, "Annotation_Barplot.pdf"), plot = anno_combined)

# ✅ TSS Distance Plot (Including Background CpGs)
tss_combined <- plotDistToTSS(
  list(
    "Hyper DMPs" = annotated_hyper_dms,
    "Hypo DMPs" = annotated_hypo_dms,
    "Hyper DMRs" = annotated_hyper_dmrs,
    "Hypo DMRs" = annotated_hypo_dmrs,
    "Background CpGs" = annotated_background
  )
)
ggsave(file.path(output_dir, "TSS_Distance_Plot.pdf"), plot = tss_combined)

cat("Updated analysis and plots generated successfully.\n")

# =========================================
# Test: Do significant DMPs drift toward 0.5 with Age?
# =========================================

# Extract beta values matrix for all CpGs
beta_mat <- as.matrix(getMeth(bs, type = "raw"))

# Keep only significant CpGs (FDR < 0.1)
sig_cpgs <- paste(DMLtest_significant$chr, DMLtest_significant$pos, sep = "_")
all_cpgs <- paste(bs@rowRanges@seqnames, bs@rowRanges@ranges@start, sep = "_")
sig_idx <- match(sig_cpgs, all_cpgs)
sig_idx <- sig_idx[!is.na(sig_idx)]
beta_sig <- beta_mat[sig_idx, , drop = FALSE]

# Compute deviation from 0.5 methylation
dev_mat <- abs(beta_sig - 0.5)

# Correlation of deviation with Age for each CpG
age_vec <- metadata$Age
cpg_corrs <- apply(dev_mat, 1, function(dev) cor(dev, age_vec, use = "complete.obs"))

# Summary statistics
mean_corr <- mean(cpg_corrs, na.rm = TRUE)
wilcox_test <- wilcox.test(cpg_corrs, mu = 0, alternative = "less") # test if correlations < 0

cat("Mean correlation between |beta-0.5| and Age across DMPs:", mean_corr, "\n")
cat("Wilcoxon signed-rank test p-value:", wilcox_test$p.value, "\n")

# =====================
# Visualization
# =====================

library(dplyr)

# Collapse CpGs to mean deviation per sample
sample_dev <- colMeans(dev_mat, na.rm = TRUE)
dev_df <- data.frame(
  Age = age_vec,
  MeanDeviation = sample_dev
)

# Plot: Mean |beta-0.5| per sample vs. Age
dev_plot <- ggplot(dev_df, aes(x = Age, y = MeanDeviation)) +
  geom_point(alpha = 0.6, color = "#2b8cbe") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  theme_minimal() +
  labs(,
    x = "Age",
    y = "Average distance from 0.5 methylation"
  )

ggsave(
  filename = file.path(output_dir, "DMP_drift_toward_mean.pdf"),
  plot = dev_plot,
  width = 7, height = 5, units = "in"
)

# Fit the linear model
lm_fit <- lm(MeanDeviation ~ Age, data = dev_df)
summary(lm_fit)
# 95% CI for slope (Age)
confint(lm_fit, "Age", level = 0.95)

# Plot: Distribution of CpG-level correlations
corr_df <- data.frame(Correlation = cpg_corrs)
corr_plot <- ggplot(corr_df, aes(x = Correlation)) +
  geom_histogram(bins = 50, fill = "#2b8cbe", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "darkred", linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Correlation between age and deviation from 0.5",
    y = "Count"
  )

ggsave(
  filename = file.path(output_dir, "DMP_correlation_distribution.pdf"),
  plot = corr_plot,
  width = 7, height = 5, units = "in"
)

# ==========================================
# Step 3b: Differential Methylation Analysis - With Cell Type Proportions
# ==========================================

# Load cell type metadata
cell_meta <- readRDS("/home/tander10/sternerlab/brain_methylation/cell_type_metadata.rds")

# Merge into existing metadata
metadata_ct <- merge(metadata, cell_meta, by = "sampid", all.x = TRUE)

# Clean up column names (spaces -> ".")
colnames(metadata_ct) <- make.names(colnames(metadata_ct))

# Ensure Sex is a factor (model.matrix will create SexM)
metadata_ct$Sex <- factor(metadata_ct$Sex, levels = c("F", "M"))

library(pheatmap)

# ================================
# Step 1: Correlation of SVs with Covariates
# ================================

# Select numeric covariates (SVs + Age + cell types)
covariates_numeric <- metadata_ct[, c("Age", "SV1", "SV2", "SV3",
                                      "astrocytes",
                                      "basket.cells",
                                      "GABAergic.neurons",
                                      "glutamatergic.neurons",
                                      "medium.spiny.neurons",
                                      "microglia",
                                      "oligodendrocytes",
                                      "vascular.cells")]

# Calculate correlation matrix
cor_matrix <- cor(covariates_numeric, use = "pairwise.complete.obs")

# Save heatmap to PDF
pdf(file.path(output_dir, "Covariate_Correlation_Heatmap.pdf"), width = 8, height = 6)
pheatmap(cor_matrix,
         main = "Correlation between SVs and Covariates",
         display_numbers = TRUE,          # show r values on heatmap
         number_format = "%.2f",
         fontsize_number = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()

cat("✅ Correlation heatmap saved to:",
    file.path(output_dir, "Covariate_Correlation_Heatmap.pdf"), "\n")

# ================================
# Step 1b: Save pairwise correlation table
# ================================

# Convert correlation matrix into long-format table
cor_table <- as.data.frame(as.table(cor_matrix))
colnames(cor_table) <- c("Var1", "Var2", "Correlation")

# Keep only upper triangle (no duplicates)
cor_table <- cor_table[as.numeric(factor(cor_table$Var1)) <
                         as.numeric(factor(cor_table$Var2)), ]

# Add flag for high correlations (|r| > 0.7)
cor_table$HighCollinearity <- ifelse(abs(cor_table$Correlation) > 0.7, "YES", "NO")

# Print top correlations sorted by |r|
cor_table_sorted <- cor_table[order(-abs(cor_table$Correlation)), ]
print(head(cor_table_sorted, 20))  # print top 20 strongest correlations

# Save full table to CSV
write.csv(cor_table_sorted,
          file = file.path(output_dir, "Covariate_Correlation_Table.csv"),
          row.names = FALSE)

cat("✅ Full correlation table saved to:",
    file.path(output_dir, "Covariate_Correlation_Table.csv"), "\n")

# --- Build design matrix manually ---
design_matrix_ct <- model.matrix(
  ~ Age + Sex + SV1 + SV2 + SV3 +
    astrocytes +
    basket.cells +
    GABAergic.neurons +
    glutamatergic.neurons +
    medium.spiny.neurons +
    microglia +
    oligodendrocytes +
    vascular.cells,
  data = metadata_ct
)
design_matrix_ct <- as.data.frame(design_matrix_ct)

# --- Run DSS ---
# IMPORTANT: model.matrix() above makes a "SexM" column, so we use SexM here
formula_ct <- ~ Age + SexM + SV1 + SV2 + SV3 +
  astrocytes +
  basket.cells +
  GABAergic.neurons +
  glutamatergic.neurons +
  medium.spiny.neurons +
  microglia +
  oligodendrocytes +
  vascular.cells

DMLfit_ct <- DMLfit.multiFactor(bs, design = design_matrix_ct, formula = formula_ct)
DMLtest_ct <- DMLtest.multiFactor(DMLfit_ct, coef = "Age")

# Save results
saveRDS(DMLtest_ct, file = file.path(output_dir, "DMLtest_results_with_SV1-3_and_CellTypes.rds"))

# Filter significant CpGs (FDR < 0.1)
DMLtest_ct_significant <- subset(DMLtest_ct, fdrs < 0.1)
summary(DMLtest_ct$fdrs < 0.1)

cat("Number of significant CpGs with cell type adjustment:", nrow(DMLtest_ct_significant), "\n")

#############################################################
# Additional Differential Methylation Models & Jackknife Tests (DSS-compatible)
#############################################################

cat("\n==============================\n")
cat(" Running Additional Models & Jackknife Tests\n")
cat("==============================\n\n")

#cell_types <- c(
#  "astrocytes", "basket.cells", "GABAergic.neurons",
#  "glutamatergic.neurons", "medium.spiny.neurons",
#  "microglia", "oligodendrocytes", "vascular.cells"
#)

cell_types <- c(
"glutamatergic.neurons", "oligodendrocytes"
)

results_summary <- data.frame(
  Model = character(),
  Removed_Covariate = character(),
  N_Significant = numeric(),
  N_Hypo = numeric(),
  Percent_Hypo = numeric(),
  stringsAsFactors = FALSE
)

# ---------------------------
# Model 1: Cell types, no SVs
# ---------------------------

cat("Running Model 1 (Cell Types Only, no SVs)...\n")

# Step 1: Construct design matrix (Sex will generate SexM)
design_ct_noSV <- model.matrix(
  as.formula(paste("~ Age + Sex +", paste(cell_types, collapse = " + "))),
  data = metadata_ct
)
design_ct_noSV <- as.data.frame(design_ct_noSV)

# Step 2: Define DSS-compatible formula (SexM now exists)
formula_ct_noSV <- as.formula(
  paste("~ Age + SexM +", paste(cell_types, collapse = " + "))
)

# Run DSS
DMLfit_ct_noSV <- DMLfit.multiFactor(bs, design = design_ct_noSV, formula = formula_ct_noSV)
DMLtest_ct_noSV <- DMLtest.multiFactor(DMLfit_ct_noSV, coef = "Age")

# Summarize results
sig_noSV <- subset(DMLtest_ct_noSV, fdrs < 0.1)
n_sig_noSV <- nrow(sig_noSV)
n_hypo_noSV <- sum(sig_noSV$stat < 0)
cat("Model 1 Results:\n")
cat("  Significant DMPs (FDR < 0.1):", n_sig_noSV, "\n")
cat("  Hypomethylated with age:", n_hypo_noSV, " (",
    round(n_hypo_noSV / max(1, n_sig_noSV) * 100, 1), "%)\n\n")

results_summary <- rbind(
  results_summary,
  data.frame(
    Model = "No_SVs",
    Removed_Covariate = "None",
    N_Significant = n_sig_noSV,
    N_Hypo = n_hypo_noSV,
    Percent_Hypo = round(n_hypo_noSV / max(1, n_sig_noSV) * 100, 1)
  )
)

# ---------------------------
# Models 2–N: Leave-one-cell-type-out ("jackknife") tests
# ---------------------------

cat("Running covariate jackknife tests (each cell type left out)...\n")

for (ct in cell_types) {
  cat("  → Running model without", ct, "...\n")
  
  vars_included <- setdiff(cell_types, ct)
  
  # Step 1: Build design matrix with Sex (creates SexM internally)
  design_jk <- model.matrix(
    as.formula(paste("~ Age + Sex + SV1 + SV2 + SV3 +", paste(vars_included, collapse = " + "))),
    data = metadata_ct
  )
  design_jk <- as.data.frame(design_jk)
  
  # Step 2: Build DSS-compatible formula (use SexM)
  #formula_jk <- as.formula(
  #  paste("~ Age + SexM + SV1 + SV2 + SV3 +", paste(vars_included, collapse = " + "))
  #)
  
  formula_jk <- as.formula(
    paste("~ Age + SexM +", paste(vars_included, collapse = " + "))
  )
  
  # Run DSS
  DMLfit_jk <- DMLfit.multiFactor(bs, design = design_jk, formula = formula_jk)
  DMLtest_jk <- DMLtest.multiFactor(DMLfit_jk, coef = "Age")
  
  # Summarize results
  sig_jk <- subset(DMLtest_jk, fdrs < 0.1)
  n_sig <- nrow(sig_jk)
  n_hypo <- sum(sig_jk$stat < 0)
  pct_hypo <- round(n_hypo / max(1, n_sig) * 100, 1)
  
  # Append to results summary
  results_summary <- rbind(
    results_summary,
    data.frame(
      Model = paste0("Jackknife_", ct),
      Removed_Covariate = ct,
      N_Significant = n_sig,
      N_Hypo = n_hypo,
      Percent_Hypo = pct_hypo
    )
  )
}

# ---------------------------
# Print summary table
# ---------------------------
cat("\n==============================\n")
cat(" Summary of Additional Models\n")
cat("==============================\n\n")
print(results_summary)

# Save to CSV
write.csv(
  results_summary,
  file = file.path(output_dir, "DMP_Model_Summary_Jackknife.csv"),
  row.names = FALSE
)

cat("\n✅ Summary table saved to DMP_Model_Summary_Jackknife.csv\n")

# ==========================================
# Step 3c: Save Significant CpGs as BED Files and Annotate
# ==========================================

# Split significant CpGs into hyper and hypo
hyper_ct <- subset(DMLtest_ct_significant, stat > 0)
hypo_ct  <- subset(DMLtest_ct_significant, stat < 0)

# Add "chr" prefix to chromosome column
hyper_ct$chr <- paste0("chr", hyper_ct$chr)
hypo_ct$chr  <- paste0("chr", hypo_ct$chr)

# Save hypermethylated CpGs as BED
hyper_bed <- data.frame(
  chrom = hyper_ct$chr,
  start = hyper_ct$pos,
  end   = hyper_ct$pos,
  name  = paste0("CpG_", seq_len(nrow(hyper_ct))),
  score = -log10(hyper_ct$fdrs),
  strand = "+"
)
write.table(hyper_bed,
            file = file.path(output_dir, "Hypermethylated_CpGs_CellTypeAdj.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Save hypomethylated CpGs as BED
hypo_bed <- data.frame(
  chrom = hypo_ct$chr,
  start = hypo_ct$pos,
  end   = hypo_ct$pos,
  name  = paste0("CpG_", seq_len(nrow(hypo_ct))),
  score = -log10(hypo_ct$fdrs),
  strand = "+"
)
write.table(hypo_bed,
            file = file.path(output_dir, "Hypomethylated_CpGs_CellTypeAdj.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("✅ Saved BED files for hyper and hypo CpGs after cell-type adjustment.\n")

# ==========================================
# Step 4: Annotation with ChIPseeker (Cell-Type Adjusted CpGs + Background)
# ==========================================

# --- Harmonize seqlevels between GRanges and TxDb ---
library(GenomeInfoDb)

# Create GRanges for background CpGs (all tested in cell-type adjusted model)
background_ct_gr <- GRanges(
  seqnames = DMLtest_ct$chr,
  ranges   = IRanges(start = DMLtest_ct$pos, end = DMLtest_ct$pos)
)

# Convert significant CpGs to GRanges
hyper_ct_gr <- GRanges(seqnames = hyper_ct$chr,
                       ranges   = IRanges(start = hyper_ct$pos, end = hyper_ct$pos))
hypo_ct_gr  <- GRanges(seqnames = hypo_ct$chr,
                       ranges   = IRanges(start = hypo_ct$pos, end = hypo_ct$pos))

# ✅ Ensure seqlevelsStyle matches TxDb (fixes "no sequence levels in common" error)
seqlevelsStyle(hyper_ct_gr)     <- seqlevelsStyle(txdb)
seqlevelsStyle(hypo_ct_gr)      <- seqlevelsStyle(txdb)
seqlevelsStyle(background_ct_gr) <- seqlevelsStyle(txdb)

# --- Run ChIPseeker annotation ---
annotated_hyper_ct <- annotatePeak(
  hyper_ct_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db"
)

annotated_hypo_ct <- annotatePeak(
  hypo_ct_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db"
)

annotated_background_ct <- annotatePeak(
  background_ct_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db"
)

# Save annotated results
saveRDS(annotated_hyper_ct,  file = file.path(output_dir, "Annotated_Hyper_CpGs_CellTypeAdj.rds"))
saveRDS(annotated_hypo_ct,   file = file.path(output_dir, "Annotated_Hypo_CpGs_CellTypeAdj.rds"))
saveRDS(annotated_background_ct, file = file.path(output_dir, "Annotated_Background_CpGs_CellTypeAdj.rds"))

# --- Plot annotation barplot including background
anno_ct_combined <- plotAnnoBar(list(
  "Hyper DMPs (CellTypeAdj)" = annotated_hyper_ct,
  "Hypo DMPs (CellTypeAdj)"  = annotated_hypo_ct,
  "Background CpGs"          = annotated_background_ct
))
ggsave(file.path(output_dir, "Annotation_Barplot_CellTypeAdj.pdf"), plot = anno_ct_combined)

# --- Plot distance to TSS including background
tss_ct_combined <- plotDistToTSS(list(
  "Hyper DMPs (CellTypeAdj)" = annotated_hyper_ct,
  "Hypo DMPs (CellTypeAdj)"  = annotated_hypo_ct,
  "Background CpGs"          = annotated_background_ct
))
ggsave(file.path(output_dir, "TSS_Distance_Plot_CellTypeAdj.pdf"), plot = tss_ct_combined)

saveRDS(metadata, file = file.path(output_dir, "metadata_filtered_86samples_with_SV1-3.rds"))
head(metadata)

cat("✅ Completed ChIPseeker annotation with background for cell-type adjusted CpGs.\n")
