# Load Required Libraries
library(Biostrings)
library(bsseq)
library(ggplot2)
library(GenomicRanges)
library(AnnotationHub)
library(ChIPseeker)
library(GenomicFeatures)
library(org.Mmu.eg.db)

# Load Required Data
mmul10_fa <- readRDS(file = "/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")
TestMeta <- read.csv(file = "/home/tander10/sternerlab/brain_methylation/Rapa_meta.csv")
bismarkBSseq_test <- readRDS(file = "/home/tander10/sternerlab/brain_methylation/bismarkBSseq_test_rapa.rds")

# Filter chromosomes to chr1-chrX
names(mmul10_fa)[1:21] -> chr
mmul10_fa <- mmul10_fa[chr]

# Add metadata to the bsseq object
pData(bismarkBSseq_test) <- TestMeta

# Filter out low coverage loci (<10x in at least 50% of samples)
coverage <- getCoverage(bismarkBSseq_test)  # Retrieve coverage matrix
min_coverage <- 10
min_samples <- ceiling(0.5 * ncol(coverage))  # At least 50% of samples

# Identify loci meeting the coverage threshold
sufficient_coverage <- rowSums(coverage >= min_coverage) >= min_samples

# Subset bsseq object to retain only high-coverage loci
bismarkBSseq_test_filtered <- bismarkBSseq_test[sufficient_coverage, ]

# Save the filtered bsseq object
saveRDS(bismarkBSseq_test_filtered, file = "/home/tander10/sternerlab/brain_methylation/Filtered_BSseq_test_rapa.rds")

# Save Percent Methylation Matrix
perc_meth <- getMeth(bismarkBSseq_test_filtered, type = "raw")
saveRDS(perc_meth, file = "/home/tander10/sternerlab/brain_methylation/Filtered_perc_meth_matrix.rds")

# PCA Visualization for Metadata
pca_result <- prcomp(t(perc_meth), center = TRUE, scale. = TRUE)
pca_df <- data.frame(pca_result$x, Group = TestMeta$Group, Year = TestMeta$Year)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = as.factor(Year))) +
  geom_point(size = 3) +
  labs(title = "PCA of Methylation Data", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("/home/tander10/sternerlab/brain_methylation/PCA_plot.pdf")

# Smooth Methylation Data
bs_smoothed <- BSmooth(bsseq = bismarkBSseq_test_filtered, verbose = TRUE)

# Differential Methylation Analysis: DMSs
tstat_DMS <- BSmooth.tstat(
  bsseq = bs_smoothed,
  group1 = TestMeta$Group == "R",
  group2 = TestMeta$Group == "C",
  estimate.var = "group",
  covariates = data.frame(Year = TestMeta$Year),
  verbose = TRUE
)

# Extract Significant CpG Sites
DMSs <- getMeth(tstat_DMS, type = "single", threshold = 4.6)
saveRDS(DMSs, file = "/home/tander10/sternerlab/brain_methylation/DMS_results.rds")

# Differential Methylation Analysis: DMRs
tstat <- BSmooth.tstat(
  bsseq = bs_smoothed,
  group1 = TestMeta$Group == "R",
  group2 = TestMeta$Group == "C",
  estimate.var = "group",
  covariates = data.frame(Year = TestMeta$Year),
  verbose = TRUE
)

dmrs <- dmrFinder(tstat, cutoff = 4.6)
saveRDS(dmrs, file = "/home/tander10/sternerlab/brain_methylation/DMR_results.rds")

# Annotate DMSs and DMRs Using ChIPseeker
# Convert DMSs and DMRs to GRanges
dms_gr <- GRanges(
  seqnames = rownames(DMSs),
  ranges = IRanges(start = DMSs$start, end = DMSs$end)
)

dmr_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges = IRanges(start = dmrs$start, end = dmrs$end)
)

# Create TxDb Object for Annotation
txdb <- makeTxDbFromEnsembl(
  organism = "Macaca mulatta",
  release = NA,
  server = "ensembldb.ensembl.org"
)

# Annotate Using ChIPseeker
annotated_dms <- annotatePeak(
  dms_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db",
  verbose = FALSE
)

annotated_dmrs <- annotatePeak(
  dmr_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db",
  verbose = FALSE
)

# Save Annotation Results
saveRDS(annotated_dms, file = "/home/tander10/sternerlab/brain_methylation/Annotated_DMS_results.rds")
saveRDS(annotated_dmrs, file = "/home/tander10/sternerlab/brain_methylation/Annotated_DMR_results.rds")

# Generate Annotation Figures
plotAnnoBar(annotated_dms) + ggtitle("DMS Annotation")
ggsave("/home/tander10/sternerlab/brain_methylation/DMS_Annotation_Barplot.pdf")

plotAnnoBar(annotated_dmrs) + ggtitle("DMR Annotation")
ggsave("/home/tander10/sternerlab/brain_methylation/DMR_Annotation_Barplot.pdf")

plotDistToTSS(annotated_dms, title = "DMS Distance to TSS")
ggsave("/home/tander10/sternerlab/brain_methylation/DMS_Distance_to_TSS.pdf")

plotDistToTSS(annotated_dmrs, title = "DMR Distance to TSS")
ggsave("/home/tander10/sternerlab/brain_methylation/DMR_Distance_to_TSS.pdf")