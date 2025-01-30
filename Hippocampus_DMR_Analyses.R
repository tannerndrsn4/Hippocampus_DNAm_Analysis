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

# Load Metadata
metadata <- read.csv(file="/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata.csv")
metadata <- as.data.frame(metadata)  # Ensure metadata is a data.frame

# Verify column names
if (!all(c("Age", "Sex") %in% colnames(metadata))) {
  stop("Metadata must contain 'Age' and 'Sex' columns.")
}

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

# Save the BSseq object in HDF5 format
library(HDF5Array)
saveHDF5SummarizedExperiment(
  bismarkBSseq,
  dir = file.path(output_dir, "bismarkBSseq_Hippocampus"),
  replace = TRUE
)

# Step 2: Filtering
bs <- filterCpGs(
  bismarkBSseq, cov = 5, perSample = 0.6,
  file = file.path(output_dir, "Filtered_BSseq_Hippocampus_5X.rds")
)

# Step 3: Differential Methylation Analysis - Site-Based
design_matrix <- model.matrix(~ Age + Sex, data = metadata)
design_matrix <- as.data.frame(design_matrix)

DMLfit <- DMLfit.multiFactor(
  bs,
  design = design_matrix,
  formula = ~ Age + SexM
)
DMLtest <- DMLtest.multiFactor(DMLfit, coef = "Age")
saveRDS(DMLtest, file = file.path(output_dir, "DMLtest_results.rds"))

# Split DMSs into hyper- and hypomethylated groups
hyper_dms <- subset(DMLtest, stat > 0 & fdrs < 0.1)
hypo_dms <- subset(DMLtest, stat < 0 & fdrs < 0.1)

# Create GRanges objects for hyper- and hypomethylated DMSs
hyper_dms_gr <- GRanges(seqnames = hyper_dms$chr, ranges = IRanges(hyper_dms$pos, hyper_dms$pos))
hypo_dms_gr <- GRanges(seqnames = hypo_dms$chr, ranges = IRanges(hypo_dms$pos, hypo_dms$pos))

# Step 4: Differential Methylation Analysis - Region-Based (Dynamic Windows)
DMRtest <- callDMR(
  DMLtest,
  p.threshold = 0.05,
  minCG = 3,
  dis.merge = 300
)
saveRDS(DMRtest, file = file.path(output_dir, "Dynamic_DMRtest_results.rds"))

# Split DMRs into hyper- and hypomethylated groups using areaStat
hyper_dmrs <- subset(DMRtest, areaStat > 0)
hypo_dmrs <- subset(DMRtest, areaStat < 0)

# Create GRanges objects for hyper- and hypomethylated DMRs
hyper_dmrs_gr <- GRanges(seqnames = hyper_dmrs$chr, ranges = IRanges(hyper_dmrs$start, hyper_dmrs$end))
hypo_dmrs_gr <- GRanges(seqnames = hypo_dmrs$chr, ranges = IRanges(hypo_dmrs$start, hypo_dmrs$end))

# Step 5: Annotation with ChIPseeker (including custom background)

# Create the txdb object
txdb <- makeTxDbFromEnsembl(
  organism = "Macaca mulatta",
  release = NA,
  server = "ensembldb.ensembl.org"
)

# Define the background: All CpG sites analyzed
background_gr <- GRanges(
  seqnames = seqnames(bs),
  ranges = IRanges(start = start(bs), end = start(bs))
)

# Annotate groups
annotated_hyper_dms <- annotatePeak(hyper_dms_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_hypo_dms <- annotatePeak(hypo_dms_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_hyper_dmrs <- annotatePeak(hyper_dmrs_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_hypo_dmrs <- annotatePeak(hypo_dmrs_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")
annotated_background <- annotatePeak(background_gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Mmu.eg.db")

# Step 6: Enhanced Visualization

# Volcano Plot for Sites
volcano_plot_sites <- ggplot(DMLtest, aes(x = stat, y = -log10(pvals))) +
  geom_point(aes(color = pvals < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red"), labels = c("Not Significant", "Significant")) +
  theme_minimal() +
  labs(title = "Volcano Plot of CpG Sites", x = "Test Statistic", y = "-log10(P-value)") +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "Volcano_Plot_Sites.pdf"), plot = volcano_plot_sites)

# Filter for significant sites
DMLtest_significant <- subset(DMLtest, fdrs < 0.1)

# Manhattan Plot for Significant Sites
manhattan_sites <- ggplot(DMLtest_significant, aes(x = pos, y = -log10(fdrs), color = factor(chr))) +
  geom_point(alpha = 0.6) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +
  theme_minimal() +
  labs(title = "Manhattan Plot of Significant CpG Sites", x = "Genomic Position", y = "-log10(FDR)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(output_dir, "Manhattan_Plot_Sites.pdf"), plot = manhattan_sites)

# Manhattan Plot for DMRs
manhattan_plot_dmrs <- ggplot(DMRtest, aes(x = start, y = abs(areaStat), color = factor(chr))) +
  geom_point(alpha = 0.6) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +
  theme_minimal() +
  labs(title = "Manhattan Plot of DMRs", x = "Genomic Position", y = "Effect Size (|AreaStat|)") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank() # Remove x-axis ticks
  )
ggsave(file.path(output_dir, "Manhattan_Plot_DMRs.pdf"), plot = manhattan_plot_dmrs)

# Annotation and Distance-to-TSS Plots
anno_combined <- plotAnnoBar(
  list(
    "Hyper DMSs" = annotated_hyper_dms,
    "Hypo DMSs" = annotated_hypo_dms,
    "Hyper DMRs" = annotated_hyper_dmrs,
    "Hypo DMRs" = annotated_hypo_dmrs,
    "Background" = annotated_background
  )
)
ggsave(file.path(output_dir, "Combined_Annotation_Barplot.pdf"), plot = anno_combined)

tss_combined <- plotDistToTSS(
  list(
    "Hyper DMSs" = annotated_hyper_dms,
    "Hypo DMSs" = annotated_hypo_dms,
    "Hyper DMRs" = annotated_hyper_dmrs,
    "Hypo DMRs" = annotated_hypo_dmrs,
    "Background" = annotated_background
  )
)
ggsave(file.path(output_dir, "Combined_Distance_to_TSS.pdf"), plot = tss_combined)

cat("All enhanced figures generated successfully.\n")
