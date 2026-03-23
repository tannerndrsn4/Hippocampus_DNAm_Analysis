# Load Required Libraries
library(ChIPseeker)
library(GenomicFeatures)
library(org.Mmu.eg.db)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(bigstatsr)
library(reshape2)
library(ggplot2)
library(bsseq)
library(comethyl)
library(txdbmaker)

# Load Metadata
metadata <- read.csv(file="/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata.csv")
metadata <- as.data.frame(metadata)  # Ensure metadata is a data.frame

# Define Paths
output_dir <- "/home/tander10/sternerlab/brain_methylation/"
dir.create(output_dir, showWarnings = FALSE)

# Verify column names
if (!all(c("Age", "Sex") %in% colnames(metadata))) {
  stop("Metadata must contain 'Age' and 'Sex' columns.")
}

# Step 1: Data Loading and Preprocessing
mmul10_fa <- readRDS(file = "/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")
names(mmul10_fa)[1:21] -> chr
mmul10_fa <- mmul10_fa[chr]
loci <- findLoci("CG", subject = mmul10_fa)

#bismarkBSseq <- read.bismark(
#  files = list.files("/home/tander10/sternerlab/brain_methylation/Bismark_output",
#                     full.names = TRUE, pattern = ".cov.gz"),
#  strandCollapse = FALSE, verbose = TRUE, BACKEND = "HDF5Array", loci = loci,
#  rmZeroCov = FALSE, dir = "/projects/sternerlab/shared/Coverage_files/hdf5",
#  replace = TRUE, BPPARAM = BiocParallel::SerialParam()
#)

# Step 2: Filtering
#bs <- filterCpGs(
#  bismarkBSseq, cov = 10, perSample = 0.6,
#  file = file.path(output_dir, "Filtered_BSseq_Hippocampus_5X.rds")
#)

#bs_file <- "/home/tander10/sternerlab/brain_methylation/results/Filtered_BSseq_Hippocampus_5X.rds"
bs_file <- "/home/tander10/sternerlab/brain_methylation/results/bs_filtered_86samples_hdf5"
library(HDF5Array)
#bs <- readRDS(bs_file)
bs <- loadHDF5SummarizedExperiment(bs_file)

# Load ARIMA Results
load("/home/tander10/sternerlab/brain_methylation/Hippocampus_arima_Region_level_SVAdj.RData")
aging_arima <- aging_arima_region

# Extract fitted and predicted values
Covfitted_nonzero <- aging_arima$fitted.nonzero
Covpredicted_nonzero <- aging_arima$predicted.nonzero
print(head(rownames(Covpredicted_nonzero)))  # Check rownames

# Format CpG indices
cpg_indices <- rownames(Covpredicted_nonzero) %>%
  # Use custom parsing
  sapply(function(name) {
    # If it starts with X followed by a number, strip the X
    name <- sub("^X(?=\\d)", "", name, perl = TRUE)
    # If it starts with _, assume it's chrX
    if (startsWith(name, "_")) {
      name <- paste0("X", name)
    }
    return(name)
  }) %>%
  as.character() %>%
  as.data.frame() %>%
  rename(index = ".") %>%
  separate(index, into = c("chrom", "start", "stop"), sep = "_", convert = TRUE) %>%
  mutate(index = paste0(chrom, "_", start, "_", stop))

# Debugging: Check index format
print("Checking cpg_indices index:")
print(head(cpg_indices$index))

# Load necessary libraries
library(dplyr)
library(tidyverse)
library(bigstatsr)
library(reshape2)
library(ggplot2)
library(cluster)  # for clusGap
library(factoextra)  # optional: for gap stat plot

# Convert matrix to Filebacked Big Matrix for chunked computation
Covpredicted_fbm <- as_FBM(t(Covpredicted_nonzero))

# Compute correlation in chunks (optimized for large matrices)
cor_blocks <- big_cor(Covpredicted_fbm, block.size = 1000)

# Convert the file-backed matrix into a numeric matrix
cor_matrix <- as.matrix(cor_blocks[])

# Convert correlation to distance matrix
clus.distance.arima <- as.dist(1 - cor_matrix)

# Perform hierarchical clustering
h.clust.arima <- hclust(clus.distance.arima, method = 'complete')

# Perform PCA to reduce dimensions before gap statistic
cat("Starting PCA reduction for gap statistic...\n")
pca_res <- prcomp(t(Covpredicted_nonzero), center = TRUE, scale. = TRUE)
# Use top 50 PCs (adjust as needed)
pca_data <- pca_res$x[, 1:50]

cat("Running gap statistic with k-means on PCA-reduced data...\n")
set.seed(123)
gap_stat <- clusGap(pca_data, FUN = kmeans, K.max = 4, B = 50)
gap_stat

# Determine optimal k
optimal_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
cat("Optimal number of clusters determined by gap statistic (PCA+kmeans):", optimal_k, "\n")

# Optional: plot and save gap statistic
pdf("/home/tander10/sternerlab/brain_methylation/Gap_statistic_plot_REGION.pdf", width = 8, height = 6)
fviz_gap_stat(gap_stat)
dev.off()

# ====================================
# 1. Gap statistic plot
# ====================================
gap_df <- as.data.frame(gap_stat$Tab)
gap_df$k <- 1:nrow(gap_df)

gap_plot <- ggplot(gap_df, aes(x = k, y = gap)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2) +
  labs(
    title = "Gap Statistic for Optimal k",
    x = "Number of clusters (k)",
    y = "Gap statistic (± SE)"
  ) +
  theme_minimal()

ggsave("GapStatisticPlot.pdf", gap_plot, width = 6, height = 4)

# ====================================
# 2. Silhouette analysis
# ====================================
# Here, use your PCA-reduced data (same as you used for gap statistic).
# Replace `pca_data` with your actual input matrix for clustering.
sil_width <- sapply(2:10, function(k) {
  km <- kmeans(pca_data, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(pca_data))
  mean(ss[, 3])  # average silhouette width
})

sil_df <- data.frame(
  k = 2:10,
  avg_silhouette = sil_width
)

sil_plot <- ggplot(sil_df, aes(x = k, y = avg_silhouette)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Average Silhouette Width Across k",
    x = "Number of clusters (k)",
    y = "Average silhouette width"
  ) +
  theme_minimal()

ggsave("SilhouettePlot.pdf", sil_plot, width = 6, height = 4)

# Print suggested k from silhouette
best_k_sil <- sil_df$k[which.max(sil_df$avg_silhouette)]
cat("Optimal number of clusters by silhouette width:", best_k_sil, "\n")

# ====================================
# 3. Gaussian Mixture Model (BIC selection)
# ====================================
library(mclust)

# Fit mixture models to your PCA-reduced data
# NOTE: If your data is very high-dimensional, it's good you already used PCA first
gmm_fit <- Mclust(pca_data, G = 1:10)  # G = range of cluster numbers

# Plot BIC for each number of clusters
pdf("GMM_BICplot.pdf", width = 6, height = 4)
plot(gmm_fit, what = "BIC")
dev.off()

# Extract the best number of clusters by BIC
best_k_gmm <- gmm_fit$G
cat("Optimal number of clusters by Gaussian Mixture Model (BIC):", best_k_gmm, "\n")

# Save the clustering result
saveRDS(h.clust.arima, file = "/home/tander10/sternerlab/brain_methylation/Cluster_dendrogram_REGION.rds")

# Cut dendrogram using optimal k
cpg_names <- rownames(Covpredicted_nonzero)
cov_gene_cluster <- tibble(
  gene = cpg_names,
  cluster = cutree(h.clust.arima, k = optimal_k)
)
cov_gene_cluster <- cov_gene_cluster %>%
  mutate(gene = gsub("^X", "", gene))

# Print cluster sizes
head(cov_gene_cluster)
print(table(cov_gene_cluster$cluster))

# Save dendrogram
pdf("/home/tander10/sternerlab/brain_methylation/Cluster_dendrogram_REGION.pdf", width = 10, height = 8)
plot(h.clust.arima, main = "Hierarchical Clustering Dendrogram (ARIMA)")
dev.off()

# Convert fitted data to dataframe
fitted_nonzero <- as.data.frame(Covfitted_nonzero)
fitted_nonzero$cluster <- cov_gene_cluster$cluster

# Split and process each cluster
for (i in sort(unique(cov_gene_cluster$cluster))) {
  cluster_names <- subset(cov_gene_cluster, cluster == i)
  cluster_data <- fitted_nonzero[fitted_nonzero$cluster == i, -ncol(fitted_nonzero)]
  rownames(cluster_data) <- cluster_names$gene
  cluster_data <- as.matrix(cluster_data)
  
  # Save CSV
  write_csv(cluster_names, paste0("/home/tander10/sternerlab/brain_methylation/cluster", i, "Regions.csv"))
  
  # Heatmap
  pdf(paste0("/home/tander10/sternerlab/brain_methylation/cluster", i, "_heatmap_REGION.pdf"))
  heatmap(cluster_data, Colv = NA, Rowv = NA)
  dev.off()
  
  # Prepare data for plotting
  cluster_df <- as.data.frame(cluster_data)
  cluster_df$Gene <- rownames(cluster_data)
  
  gene_expression_melted <- melt(cluster_df, id.vars = "Gene", variable.name = "Age", value.name = "Expression")
  gene_expression_melted$Age <- as.numeric(gsub("X", "", gene_expression_melted$Age))
  
  summary_stats <- gene_expression_melted %>%
    group_by(Age) %>%
    summarise(mean_expr = mean(Expression, na.rm = TRUE),
              sd_expr = sd(Expression, na.rm = TRUE),
              .groups = 'drop')
  
  # Store mins and maxs for each cluster
  smoothed_y_mins <- c()
  smoothed_y_maxs <- c()

  # First pass: compute smoothed y range for each cluster
  for (i in sort(unique(cov_gene_cluster$cluster))) {
    cluster_names <- subset(cov_gene_cluster, cluster == i)
    cluster_data <- fitted_nonzero[fitted_nonzero$cluster == i, -ncol(fitted_nonzero)]
    rownames(cluster_data) <- cluster_names$gene
    cluster_data <- as.matrix(cluster_data)
    
    cluster_df <- as.data.frame(cluster_data)
    cluster_df$Gene <- rownames(cluster_data)
    
    gene_expression_melted <- melt(cluster_df, id.vars = "Gene", variable.name = "Age", value.name = "Expression")
    gene_expression_melted$Age <- as.numeric(gsub("X", "", gene_expression_melted$Age))
    
    summary_stats <- gene_expression_melted %>%
      group_by(Age) %>%
      summarise(mean_expr = mean(Expression, na.rm = TRUE),
                .groups = 'drop')
    
    # Fit LOESS and predict
    loess_fit <- loess(mean_expr ~ Age, data = summary_stats, span = 0.75)
    smoothed_vals <- predict(loess_fit, newdata = summary_stats$Age)
    
    smoothed_y_mins <- c(smoothed_y_mins, min(smoothed_vals, na.rm = TRUE))
    smoothed_y_maxs <- c(smoothed_y_maxs, max(smoothed_vals, na.rm = TRUE))
  }

  library(pracma)  # For numerical differentiation

  # Compute global y-limits with buffer
  global_y_min <- min(smoothed_y_mins) - 0.005
  global_y_max <- max(smoothed_y_maxs) + 0.005

  # Now generate the plots with the shared y-axis range
  for (i in sort(unique(cov_gene_cluster$cluster))) {
    cluster_names <- subset(cov_gene_cluster, cluster == i)
    cluster_data <- fitted_nonzero[fitted_nonzero$cluster == i, -ncol(fitted_nonzero)]
    rownames(cluster_data) <- cluster_names$gene
    cluster_data <- as.matrix(cluster_data)
    
    cluster_df <- as.data.frame(cluster_data)
    cluster_df$Gene <- rownames(cluster_data)
    
    gene_expression_melted <- melt(cluster_df, id.vars = "Gene", variable.name = "Age", value.name = "Expression")
    gene_expression_melted$Age <- as.numeric(gsub("X", "", gene_expression_melted$Age))
    
    summary_stats <- gene_expression_melted %>%
      group_by(Age) %>%
      summarise(mean_expr = mean(Expression, na.rm = TRUE),
                sd_expr = sd(Expression, na.rm = TRUE),
                .groups = 'drop')

    # Fit LOESS model
    loess_fit <- loess(mean_expr ~ Age, data = summary_stats, span = 0.75)

    # Predict on fine-grained age grid
    age_grid <- seq(min(summary_stats$Age), max(summary_stats$Age), by = 0.1)
    predicted_vals <- predict(loess_fit, newdata = data.frame(Age = age_grid))

    # Identify max and min expression points
    max_idx <- which.max(predicted_vals)
    min_idx <- which.min(predicted_vals)

    max_age <- age_grid[max_idx]
    min_age <- age_grid[min_idx]
    max_val <- predicted_vals[max_idx]
    min_val <- predicted_vals[min_idx]

    cat(paste0("Cluster ", i, ":\n"))
    cat(paste0("  Maximum methylation at age ", round(max_age, 2), " (", round(max_val, 3), ")\n"))
    cat(paste0("  Minimum methylation at age ", round(min_age, 2), " (", round(min_val, 3), ")\n"))

    # Create smooth plot
    avg_smooth <- ggplot(summary_stats, aes(x = Age, y = mean_expr)) +
      geom_smooth(se = FALSE, color = "#2b8cbe", method = "loess", span = 0.75, linewidth = 1.5) +
      theme_classic() +
      labs(title = paste("Methylation % Trajectory (Cluster", i, ")"),
           x = "Age",
           y = "Methylation % trajectory") +
      coord_cartesian(ylim = c(global_y_min, global_y_max)) +
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 16))

    print(avg_smooth)
    ggsave(filename = paste0("/home/tander10/sternerlab/brain_methylation/CovCluster", i, "_smooth_REGION_plot.pdf"),
           plot = avg_smooth)
  }
         
  # Smoothed plot with ribbon
  
  #avg_smooth <- ggplot(summary_stats, aes(x = Age, y = mean_expr)) +
    #geom_ribbon(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
                #fill = "gray80", alpha = 0.5) +
    #geom_smooth(se = FALSE, color = "#2b8cbe", method = "loess", size = 1.2) +
    #theme_classic() +
    #labs(title = paste("Methylation % Trajectory (Cluster", i, ")"),
         #x = "Age",
         #y = "Methylation % Trajectory") +
    #ylim(0.55, 0.85) +
    #theme(axis.title = element_text(size = 20),
          #axis.text = element_text(size = 16))
  
  #print(avg_smooth)
  #ggsave(filename = paste0("/home/tander10/sternerlab/brain_methylation/CovCluster", i, "_smooth_expression_plot.pdf"),
         #plot = avg_smooth)
}


# ---- Generate GRanges objects dynamically for each cluster ----

# Merge cov_gene_cluster (with cluster labels) and cpg_indices (with genomic positions)
cov_gene_cluster <- cov_gene_cluster %>%
  mutate(gene = sapply(gene, function(name) {
    # Strip "X" if followed by a digit (e.g., "X1_...")
    name <- sub("^X(?=\\d)", "", name, perl = TRUE)
    # If it starts with "_", it's chrX
    if (startsWith(name, "_")) {
      name <- paste0("X", name)
    }
    return(name)
  }))

head(cov_gene_cluster)
tail(cov_gene_cluster)
head(cpg_indices)
tail(cpg_indices)

cluster_annot_data <- cov_gene_cluster %>%
  left_join(cpg_indices, by = c("gene" = "index"))

summary(is.na(cluster_annot_data$start))  # TRUE means there's at least one NA
summary(is.na(cluster_annot_data$chrom))  # Useful for confirming join match
head(cluster_annot_data)

# Check for any missing joins
if (any(is.na(cluster_annot_data$chrom))) {
  warning("Some CpGs could not be matched to genomic coordinates. Check cpg_indices and cov_gene_cluster consistency.")
}

# Create list to hold GRanges objects for each cluster
cluster_gr_list <- list()

for (i in sort(unique(cluster_annot_data$cluster))) {
  this_cluster <- cluster_annot_data %>% filter(cluster == i)
  
  gr <- GRanges(
    seqnames = this_cluster$chrom,
    ranges = IRanges(start = this_cluster$start, end = this_cluster$stop)
  )
  
  cluster_gr_list[[paste0("Cluster", i)]] <- gr
}

# ---- Annotate each cluster ----

# Create txdb object
txdb <- txdbmaker::makeTxDbFromEnsembl(
  organism = "Macaca mulatta",
  release = NA,
  server = "ensembldb.ensembl.org"
)

saveDb(txdb, file.path(output_dir, "Macaca_mulatta_txdb.sqlite"))
#txdb <- loadDb("/home/tander10/sternerlab/brain_methylation/Macaca_mulatta_txdb.sqlite")

# Create list to hold annotation results
annotated_clusters <- list()

for (name in names(cluster_gr_list)) {
  annotated <- annotatePeak(
    cluster_gr_list[[name]],
    TxDb = txdb,
    tssRegion = c(-3000, 3000),
    annoDb = "org.Mmu.eg.db"
  )
  
  annotated_clusters[[name]] <- annotated
}
head(annotated_clusters)

# Save all annotated cluster objects
save(annotated_clusters, file = file.path(output_dir, "Annotated_REGION_Clusters.RData"))
cat("Annotated REGION cluster results saved successfully in", file.path(output_dir, "Annotated_REGION_Clusters.RData"), "\n")

# ---- Annotate background CpGs ----

background_cpgs <- rowRanges(bs)
background_cpgs_gr <- GRanges(seqnames = seqnames(background_cpgs), ranges = ranges(background_cpgs))

annotated_background <- annotatePeak(
  background_cpgs_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db"
)

saveRDS(annotated_background, file = file.path(output_dir, "annotated_background.rds"))
annotated_background@annoStat

# ---- Create GRanges for ALL age-associated regions ----

all_regions_gr <- GRanges(
  seqnames = cluster_annot_data$chrom,
  ranges = IRanges(start = cluster_annot_data$start,
                   end = cluster_annot_data$stop)
)

# Annotate them
annotated_all_regions <- annotatePeak(
  all_regions_gr,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mmu.eg.db"
)

# Save annotated regions object
saveRDS(annotated_all_regions,
        file = file.path(output_dir, "All_AgeAssociated_Regions_ChIPseeker_Annotation.rds"))

# ---- Plot annotations ----

extract_annotation_df <- function(peak_obj, label){

  df <- as.data.frame(peak_obj@anno)

  df$Feature <- df$annotation

  # ---- Promoter / Intergenic ----
  df$Feature[grepl("Promoter", df$Feature)] <- "Promoter"
  df$Feature[grepl("Distal Intergenic", df$Feature)] <- "Distal Intergenic"

  # ---- UTRs ----
  df$Feature[grepl("5' UTR", df$Feature)] <- "5' UTR"
  df$Feature[grepl("3' UTR", df$Feature)] <- "3' UTR"

  # ---- FIRST exon/intron (must come before general rules) ----
  df$Feature[grepl("exon 1 of", df$Feature, ignore.case = TRUE)] <- "1st Exon"
  df$Feature[grepl("intron 1 of", df$Feature, ignore.case = TRUE)] <- "1st Intron"

  # ---- Other exons/introns ----
  df$Feature[grepl("^Exon", df$Feature) & !grepl("exon 1 of", df$annotation, ignore.case = TRUE)] <- "Other Exon"
  df$Feature[grepl("^Intron", df$Feature) & !grepl("intron 1 of", df$annotation, ignore.case = TRUE)] <- "Other Intron"

  # ---- Downstream ----
  df$Feature[grepl("Downstream", df$Feature)] <- "Downstream"

  df$Dataset <- label

  return(df[,c("Feature","Dataset")])
}

anno_list <- list()

anno_list[[1]] <- extract_annotation_df(annotated_all_regions, "All age-associated regions")

for(i in seq_along(annotated_clusters)){
  anno_list[[i+1]] <- extract_annotation_df(annotated_clusters[[i]], paste0("Cluster",i))
}

anno_list[[length(anno_list)+1]] <- extract_annotation_df(annotated_background, "All background regions")

anno_df <- bind_rows(anno_list)

anno_summary <- anno_df %>%
  group_by(Dataset, Feature) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(Dataset) %>%
  mutate(Percent = 100*n/sum(n))

anno_summary$Dataset <- factor(
  anno_summary$Dataset,
  levels = c(
    "All background regions",
    "Cluster 4",
    "Cluster 3",
    "Cluster 2",
    "Cluster 1",
    "All age-associated regions"
  )
)

feature_colors <- c(
  "Distal Intergenic" = "#593787",
  "Promoter" = "#56ae6a",
  "5' UTR" = "#aaa23f",
  "1st Exon" = "#6973d8",
  "Other Exon" = "#c074cc",
  "1st Intron" = "#b85f36",
  "Other Intron" = "#bc4862",
  "3' UTR" = "#5794d7",
  "Downstream" = "#b5508f"
)

anno_plot <- ggplot(anno_summary,
                    aes(x=Dataset, y=Percent, fill=Feature)) +
  geom_bar(stat="identity", position="fill") +
  coord_flip() +
  scale_fill_manual(values=feature_colors) +
  ylab("Percent (%)") +
  xlab("") +
  theme_classic() +
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    legend.title = element_blank()
  )

ggsave(file.path(output_dir,"REGION_Annotation_Barplot.pdf"),
       anno_plot,
       width=8,
       height=5)

# ---- Export clusters and background as BED files ----

# Function to export GRanges as BED with chr prefix
export_as_bed <- function(gr, file) {
  # Make sure seqnames have "chr" prefix
  seqlevelsStyle(gr) <- "UCSC"  # ensures chr1, chr2, chrX format
  df <- as.data.frame(gr)[, c("seqnames", "start", "end")]
  write.table(df, file = file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# Export each cluster
for (name in names(cluster_gr_list)) {
  out_file <- file.path(output_dir, paste0(name, "_regions.bed"))
  export_as_bed(cluster_gr_list[[name]], out_file)
  cat("Exported:", out_file, "\n")
}

# Export background CpGs
seqlevelsStyle(background_cpgs_gr) <- "UCSC"
background_bed <- file.path(output_dir, "Background_CpGs.bed")
export_as_bed(background_cpgs_gr, background_bed)
cat("Exported:", background_bed, "\n")

library(rtracklayer)
library(GenomeInfoDb)

# ---- Define ARIMA Input Regions EXACTLY as in ARIMA analysis ----

# Load residual_region (from ARIMA pipeline)
residual_region <- readRDS("/home/tander10/sternerlab/brain_methylation/Percent_Methylated_Regions_Residuals_SexSVAdj_Limma.rds")

# Extract region IDs from rownames (format: chr_start_end)
region_ids <- rownames(residual_region)

# Split into components
region_split <- strsplit(region_ids, "_")
region_df <- do.call(rbind, region_split)
colnames(region_df) <- c("chr", "start", "end")

# Convert to data.frame
region_df <- as.data.frame(region_df, stringsAsFactors = FALSE)
region_df$start <- as.integer(region_df$start)
region_df$end   <- as.integer(region_df$end)
region_df$RegionID <- region_ids

# Build GRanges for background (mmul10 coordinates)
arima_regions_gr <- GRanges(
  seqnames = region_df$chr,
  ranges   = IRanges(start = region_df$start, end = region_df$end),
  RegionID = region_df$RegionID
)

cat("Loaded ARIMA background regions:", length(arima_regions_gr), "\n")

# ---- Liftover + Sync Foreground & Background ----

# Import chain
chain_file <- "/home/tander10/sternerlab/brain_methylation/rheMac10ToHg38.over.chain"
chain <- rtracklayer::import.chain(chain_file)

# Function to run liftover + report losses
liftover_with_report <- function(gr, chain, name) {
  original_n <- length(gr)
  lifted <- liftOver(gr, chain)
  lifted_unlisted <- unlist(lifted)
  lifted_n <- length(lifted_unlisted)
  cat(sprintf("%s: %d original regions, %d lifted regions (dropped %d)\n",
              name, original_n, lifted_n, original_n - lifted_n))
  return(lifted_unlisted)
}

# Liftover background (ARIMA regions)
seqlevelsStyle(arima_regions_gr) <- "UCSC"
lifted_background <- liftover_with_report(arima_regions_gr, chain, "Background")

# Liftover each cluster
lifted_clusters <- list()
for (name in names(cluster_gr_list)) {
  seqlevelsStyle(cluster_gr_list[[name]]) <- "UCSC"
  lifted_clusters[[name]] <- liftover_with_report(cluster_gr_list[[name]], chain, name)
}

# ---- Ensure Foreground ⊆ Background ----
# Keep only cluster regions that exist in background after liftover
bg_key <- paste0(seqnames(lifted_background), "_", start(lifted_background), "_", end(lifted_background))

for (name in names(lifted_clusters)) {
  cl_key <- paste0(seqnames(lifted_clusters[[name]]), "_", start(lifted_clusters[[name]]), "_", end(lifted_clusters[[name]]))
  keep_idx <- cl_key %in% bg_key
  dropped <- sum(!keep_idx)
  cat(sprintf("%s: %d regions dropped because not in background\n", name, dropped))
  lifted_clusters[[name]] <- lifted_clusters[[name]][keep_idx]
}

# ---- Export synced sets ----
export_as_bed(lifted_background, file.path(output_dir, "ARIMA_input_regions_hg38.bed"))

for (name in names(lifted_clusters)) {
  out_file <- file.path(output_dir, paste0(name, "_regions_hg38.bed"))
  export_as_bed(lifted_clusters[[name]], out_file)
  cat("Exported:", out_file, "\n")
}

# ====================================
# Plot smoothed methylation trajectories for specific regions
# ====================================

# --- Load libraries ---
library(ggplot2)
library(dplyr)

# --- Helper function to plot methylation vs age for a region ---
plot_region_smooth <- function(region_id, methylation_matrix) {
  # Extract and reshape region methylation
  df <- as.data.frame(t(methylation_matrix[region_id, , drop = FALSE]))
  
  # Rename column
  colnames(df) <- "Methylation"
  
  # Extract numeric age from column names (assuming colnames like "X10", "X20", etc.)
  df$Age <- as.numeric(gsub("X", "", rownames(df)))
  
  # Remove missing values
  df <- df %>% filter(!is.na(Methylation) & !is.na(Age))
  
  # Skip if too few points
  if (nrow(df) < 5 || length(unique(df$Age)) < 3) {
    cat(paste0("⚠️ Skipping ", region_id, " — insufficient data for smoothing.\n"))
    return(NULL)
  }
  
  # --- Fit LOESS model manually to calculate smoothed range ---
  loess_fit <- loess(Methylation ~ Age, data = df, span = 0.75)
  age_grid <- seq(min(df$Age), max(df$Age), by = 0.1)
  smooth_vals <- predict(loess_fit, newdata = data.frame(Age = age_grid))
  
  # Compute dynamic y-limits with ±0.2 buffer
  ymin <- min(smooth_vals, na.rm = TRUE) - 0.05
  ymax <- max(smooth_vals, na.rm = TRUE) + 0.05
  
  # --- Build clean plot ---
  p <- ggplot(df, aes(x = Age, y = Methylation)) +
    geom_smooth(method = "loess", color = "#2b8cbe", se = FALSE, linewidth = 1.3) +
    theme_classic(base_size = 14) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(
      title = paste("Smoothed Methylation Trajectory:", region_id),
      x = "Age (years)",
      y = "Percent Methylation"
    ) +
    scale_x_continuous(breaks = c(10, 20, 30)) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# --- Define target regions ---
target_regions <- c("7_77879817_77879847", "7_135272269_135272542", "1_203983149_203983234")

# --- Generate and save plots ---
plots <- list()

for (region in target_regions) {
  if (region %in% rownames(residual_region)) {
    cat(paste0("Plotting ", region, "...\n"))
    p <- plot_region_smooth(region, residual_region)
    
    if (!is.null(p)) {
      plots[[region]] <- p
      
      # Save plots
      out_pdf <- paste0("/home/tander10/sternerlab/brain_methylation/", region, "_Methylation_Trajectory.pdf")
      out_png <- paste0("/home/tander10/sternerlab/brain_methylation/", region, "_Methylation_Trajectory.png")
      
      ggsave(out_pdf, p, width = 6, height = 4)
      
      cat("✅ Saved:", region, "\n")
    }
  } else {
    cat(paste0("⚠️ Region ", region, " not found in residual_region.\n"))
  }
}


# ===============================
# SEX-SPECIFIC ARIMA CLUSTERING + TRAJECTORY PLOTS (REGION LEVEL)
# Original-style plotting (no age_map joins)
# ===============================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(reshape2)
  library(bigstatsr)
  library(cluster)     # clusGap, silhouette
  library(factoextra)  # fviz_gap_stat
  library(mclust)      # GMM
})

# ---- Inputs ----
arima_rdata_F <- "/home/tander10/sternerlab/brain_methylation/results_arima_region_sex/ARIMA_Region_F/Hippocampus_arima_Region_level_F_SVAdj.RData"
arima_rdata_M <- "/home/tander10/sternerlab/brain_methylation/results_arima_region_sex/ARIMA_Region_M/Hippocampus_arima_Region_level_M_SVAdj.RData"

output_root <- "/home/tander10/sternerlab/brain_methylation/results_arima_region_sex_clustering"
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

sex_colors <- c(F = "#b561b8", M = "#6980ce")

set.seed(123)
K_MAX_GAP <- 4
GAP_B <- 50
PCA_NPCS <- 50
GMM_MAX_G <- 10
LOESS_SPAN <- 0.75

safe_npcs <- function(x, requested = 50) {
  # x is matrix rows=regions, cols=ages
  # prcomp on t(x): rows = ages, cols = regions (or vice versa depending)
  n <- min(nrow(t(x)) - 1, ncol(t(x)) - 1, requested)
  max(2, n)
}

# ---- Core runner (one sex) ----
run_one_sex <- function(arima_rdata, sex_label, outdir, line_color_hex) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  cat("\n====================================\n")
  cat("Running clustering pipeline for sex:", sex_label, "\n")
  cat("ARIMA file:", arima_rdata, "\n")
  cat("Output dir:", outdir, "\n")
  cat("====================================\n\n")

  e <- new.env()
  load(arima_rdata, envir = e)

  # Grab ARIMA object from RData
  if (exists("aging_arima_region_sx", envir = e)) {
    arima_obj <- get("aging_arima_region_sx", envir = e)
  } else if (exists("aging_arima_region", envir = e)) {
    arima_obj <- get("aging_arima_region", envir = e)
  } else if (length(ls(envir = e)) == 1) {
    arima_obj <- get(ls(envir = e)[1], envir = e)
    warning("Using the only object found in RData: ", ls(envir = e)[1])
  } else {
    stop("Could not find ARIMA object in: ", arima_rdata)
  }

  Covfitted_nonzero    <- arima_obj$fitted.nonzero
  Covpredicted_nonzero <- arima_obj$predicted.nonzero
  rownames(Covpredicted_nonzero) <- rownames(Covfitted_nonzero)
  
  # =========================
  # DIAGNOSTICS: fitted vs predicted
  # =========================
  cat("\n--- DIAGNOSTICS:", sex_label, "---\n")

  cat("predicted.nonzero dims:", paste(dim(Covpredicted_nonzero), collapse = " x "), "\n")
  cat("fitted.nonzero    dims:", paste(dim(Covfitted_nonzero), collapse = " x "), "\n")

  # Rownames previews
  cat("\nRownames (predicted) head:\n")
  print(utils::head(rownames(Covpredicted_nonzero)))
  cat("Rownames (predicted) tail:\n")
  print(utils::tail(rownames(Covpredicted_nonzero)))

  cat("\nRownames (fitted) head:\n")
  print(utils::head(rownames(Covfitted_nonzero)))
  cat("Rownames (fitted) tail:\n")
  print(utils::tail(rownames(Covfitted_nonzero)))

  # Identity checks (set + order)
  same_row_order <- identical(rownames(Covpredicted_nonzero), rownames(Covfitted_nonzero))
  cat("\nDo predicted and fitted have IDENTICAL rownames (same order)? ", same_row_order, "\n")

  # Set comparisons
  pred_ids <- rownames(Covpredicted_nonzero)
  fit_ids  <- rownames(Covfitted_nonzero)

  shared_ids <- intersect(pred_ids, fit_ids)
  pred_only  <- setdiff(pred_ids, fit_ids)
  fit_only   <- setdiff(fit_ids, pred_ids)

  cat("Shared RegionIDs:", length(shared_ids), "\n")
  cat("Pred-only RegionIDs:", length(pred_only), "\n")
  cat("Fit-only RegionIDs:", length(fit_only), "\n")

  if (length(pred_only) > 0) {
    cat("\nExample pred-only RegionIDs (first 10):\n")
    print(utils::head(pred_only, 10))
  }
  if (length(fit_only) > 0) {
    cat("\nExample fit-only RegionIDs (first 10):\n")
    print(utils::head(fit_only, 10))
  }

  # Column name checks
  pred_cols <- colnames(Covpredicted_nonzero)
  fit_cols  <- colnames(Covfitted_nonzero)

  cat("\npredicted.nonzero colnames head:", paste(utils::head(pred_cols), collapse = ", "), "\n")
  cat("fitted.nonzero    colnames head:", paste(utils::head(fit_cols), collapse = ", "), "\n")

  same_col_order <- identical(pred_cols, fit_cols)
  cat("Do predicted and fitted have IDENTICAL colnames (same order)? ", same_col_order, "\n")

  # If colnames differ, show quick set diffs
  if (!same_col_order) {
    cat("Shared colnames:", length(intersect(pred_cols, fit_cols)), "\n")
    cat("Pred-only colnames:", length(setdiff(pred_cols, fit_cols)), "\n")
    cat("Fit-only colnames:", length(setdiff(fit_cols, pred_cols)), "\n")
  }

  cat("--- END DIAGNOSTICS ---\n\n")

  if (is.null(Covpredicted_nonzero) || nrow(Covpredicted_nonzero) == 0) {
    stop("predicted.nonzero is empty for sex ", sex_label)
  }
  if (is.null(Covfitted_nonzero) || nrow(Covfitted_nonzero) == 0) {
    stop("fitted.nonzero is empty for sex ", sex_label)
  }

  # ---- CRITICAL: force fitted/predicted to share the same RegionIDs ----
  common_regions <- intersect(rownames(Covpredicted_nonzero), rownames(Covfitted_nonzero))
  if (length(common_regions) < 10) {
    stop("Too few shared regions between predicted and fitted for sex ", sex_label,
         ". Shared=", length(common_regions))
  }
  Covpredicted_nonzero <- Covpredicted_nonzero[common_regions, , drop = FALSE]
  Covfitted_nonzero    <- Covfitted_nonzero[common_regions, , drop = FALSE]

  cat("Regions (shared):", length(common_regions), "\n")
  cat("Eval ages (pred):", ncol(Covpredicted_nonzero), "\n")
  cat("Eval ages (fit): ", ncol(Covfitted_nonzero), "\n")

  # =========================
  # 1) PCA + k diagnostics
  # =========================
  npcs <- safe_npcs(Covpredicted_nonzero, PCA_NPCS)
  cat("Running PCA (keeping", npcs, "PCs)...\n")
  pca_res <- prcomp(t(Covpredicted_nonzero), center = TRUE, scale. = TRUE)
  pca_data <- pca_res$x[, 1:npcs, drop = FALSE]

  cat("Running gap statistic (k=1..", K_MAX_GAP, ", B=", GAP_B, ")...\n", sep = "")
  gap_stat <- clusGap(pca_data, FUN = kmeans, K.max = K_MAX_GAP, B = GAP_B)
  optimal_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
  cat("Optimal k by gap statistic:", optimal_k, "\n")

  pdf(file.path(outdir, paste0("Gap_statistic_", sex_label, ".pdf")), width = 8, height = 6)
  print(fviz_gap_stat(gap_stat))
  dev.off()

  # Silhouette
  cat("Running silhouette analysis...\n")
  sil_width <- sapply(2:10, function(k) {
    km <- kmeans(pca_data, centers = k, nstart = 25)
    ss <- silhouette(km$cluster, dist(pca_data))
    mean(ss[, 3])
  })
  sil_df <- data.frame(k = 2:10, avg_silhouette = sil_width)
  best_k_sil <- sil_df$k[which.max(sil_df$avg_silhouette)]
  cat("Optimal k by silhouette:", best_k_sil, "\n")

  ggsave(
    file.path(outdir, paste0("SilhouettePlot_", sex_label, ".pdf")),
    ggplot(sil_df, aes(k, avg_silhouette)) + geom_line() + geom_point() +
      theme_minimal() + labs(title=paste0("Silhouette (", sex_label, ")")),
    width = 6, height = 4
  )

  # GMM
  cat("Running GMM (mclust) BIC selection...\n")
  gmm_fit <- Mclust(pca_data, G = 1:GMM_MAX_G)
  best_k_gmm <- gmm_fit$G
  cat("Optimal k by GMM/BIC:", best_k_gmm, "\n")

  pdf(file.path(outdir, paste0("GMM_BICplot_", sex_label, ".pdf")), width = 6, height = 4)
  plot(gmm_fit, what = "BIC")
  dev.off()

  # =========================
  # 2) Hierarchical clustering on predicted trajectories
  # =========================
  cat("Computing correlation (bigstatsr) + hierarchical clustering...\n")
  Covpredicted_fbm <- as_FBM(t(Covpredicted_nonzero))
  cor_blocks <- big_cor(Covpredicted_fbm, block.size = 1000)
  cor_matrix <- as.matrix(cor_blocks[])
  dist_mat <- as.dist(1 - cor_matrix)

  hclust_obj <- hclust(dist_mat, method = "complete")
  saveRDS(hclust_obj, file.path(outdir, paste0("Cluster_dendrogram_", sex_label, ".rds")))

  pdf(file.path(outdir, paste0("Cluster_dendrogram_", sex_label, ".pdf")), width = 10, height = 8)
  plot(hclust_obj, main = paste0("Hierarchical Clustering Dendrogram (", sex_label, ")"))
  dev.off()

  clusters <- cutree(hclust_obj, k = optimal_k)
  names(clusters) <- rownames(Covpredicted_nonzero)

  cat("Cluster sizes:\n")
  print(table(clusters))

  cluster_df <- tibble(RegionID = names(clusters), cluster = as.integer(clusters))
  write_csv(cluster_df, file.path(outdir, paste0("Region_clusters_", sex_label, "_k", optimal_k, ".csv")))

  # =========================
  # 3) Heatmaps + mean trajectory plots (ORIGINAL STYLE)
  # =========================
  fitted_mat <- as.matrix(Covfitted_nonzero)

  # Helper: make per-cluster long df, parse Age from colnames directly (like your original)
  cluster_long <- function(k) {
    regs <- names(clusters)[clusters == k]
    regs <- intersect(regs, rownames(fitted_mat))  # extra safety
    if (length(regs) < 2) return(NULL)

    mat_k <- fitted_mat[regs, , drop = FALSE]
    df <- as.data.frame(mat_k)
    df$RegionID <- rownames(mat_k)

    m <- melt(df, id.vars="RegionID", variable.name="Age", value.name="Methylation")
    m$Age <- suppressWarnings(as.numeric(sub("^X", "", as.character(m$Age))))
    m <- m %>% filter(is.finite(Age), is.finite(Methylation))
    if (nrow(m) == 0) return(NULL)
    m
  }

  # First pass: compute global y-limits from loess-smoothed mean across clusters
  smoothed_y_mins <- c()
  smoothed_y_maxs <- c()

  for (k in sort(unique(clusters))) {
    m <- cluster_long(k)
    if (is.null(m)) next

    ss <- m %>%
      group_by(Age) %>%
      summarise(mean_meth = mean(Methylation, na.rm=TRUE), .groups="drop") %>%
      arrange(Age)

    if (nrow(ss) < 5 || length(unique(ss$Age)) < 3) next

    lo <- tryCatch(loess(mean_meth ~ Age, data=ss, span=LOESS_SPAN), error=function(e) NULL)
    if (is.null(lo)) next

    age_grid <- seq(min(ss$Age), max(ss$Age), by = 0.1)
    pred <- predict(lo, newdata=data.frame(Age=age_grid))
    if (all(is.na(pred))) next

    smoothed_y_mins <- c(smoothed_y_mins, min(pred, na.rm=TRUE))
    smoothed_y_maxs <- c(smoothed_y_maxs, max(pred, na.rm=TRUE))
  }

  if (length(smoothed_y_mins) == 0) {
    warning("No valid loess fits for global y-limits; using fitted_mat range.")
    global_y_min <- min(fitted_mat, na.rm=TRUE) - 0.005
    global_y_max <- max(fitted_mat, na.rm=TRUE) + 0.005
  } else {
    global_y_min <- min(smoothed_y_mins, na.rm=TRUE) - 0.005
    global_y_max <- max(smoothed_y_maxs, na.rm=TRUE) + 0.005
  }

  n_saved <- 0

  for (k in sort(unique(clusters))) {

    # Heatmap
    regs <- names(clusters)[clusters == k]
    regs <- intersect(regs, rownames(fitted_mat))
    if (length(regs) >= 2) {
      pdf(file.path(outdir, paste0("Cluster", k, "_heatmap_", sex_label, ".pdf")))
      heatmap(fitted_mat[regs, , drop=FALSE], Colv=NA, Rowv=NA, scale="none")
      dev.off()
    }

    # Mean trajectory plot
    m <- cluster_long(k)
    if (is.null(m)) next

    ss <- m %>%
      group_by(Age) %>%
      summarise(mean_meth = mean(Methylation, na.rm=TRUE),
                sd_meth   = sd(Methylation, na.rm=TRUE),
                .groups="drop") %>%
      arrange(Age)

    if (nrow(ss) < 5 || length(unique(ss$Age)) < 3) next

    p <- ggplot(ss, aes(x=Age, y=mean_meth)) +
      geom_smooth(se=FALSE, method="loess", span=LOESS_SPAN,
                  linewidth=1.5, color=line_color_hex) +
      theme_classic() +
      coord_cartesian(ylim=c(global_y_min, global_y_max)) +
      labs(title=paste0("Mean methylation trajectory (", sex_label, ") — Cluster ", k),
           x="Age", y="Mean methylation (fitted)")

    outfile <- file.path(outdir, paste0("Cluster", k, "_meanTrajectory_", sex_label, ".pdf"))
    ggsave(outfile, p, width=6.5, height=4.5)
    cat("✅ Saved:", outfile, "\n")
    n_saved <- n_saved + 1
  }

  if (n_saved == 0) {
    stop("No trajectory PDFs were saved for sex ", sex_label,
         ". This usually means Age parsing failed or summary_stats had too few unique ages.")
  }

  # Save k summary
  k_summary <- tibble(
    sex = sex_label,
    k_gap = optimal_k,
    k_silhouette = best_k_sil,
    k_gmm = best_k_gmm,
    k_used = optimal_k,
    n_regions = nrow(Covpredicted_nonzero)
  )
  write_csv(k_summary, file.path(outdir, paste0("k_summary_", sex_label, ".csv")))

  cat("\n✅ Done for sex:", sex_label, "\n")

  invisible(list(
    cluster_df = cluster_df,
    clusters = clusters,
    region_ids = names(clusters)
  ))
}

# ---- Run ----
res_F <- run_one_sex(arima_rdata_F, "F", file.path(output_root, "Female"), sex_colors["F"])
res_M <- run_one_sex(arima_rdata_M, "M", file.path(output_root, "Male"),   sex_colors["M"])

# ===============================
# Overlap comparisons (Sex vs Sex, Sex vs Original)
# ===============================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(readr)
  library(dplyr)
})

regionid_to_granges <- function(region_ids, add_chr = TRUE) {
  region_ids <- as.character(region_ids)
  parts <- do.call(rbind, strsplit(region_ids, "_"))
  if (ncol(parts) != 3) stop("RegionIDs must look like 'chr_start_end' or '1_start_end'.")
  df <- as.data.frame(parts, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end")
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  if (add_chr && !grepl("^chr", df$chr[1])) df$chr <- paste0("chr", df$chr)
  GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
}

overlap_percent_any <- function(A, B) {
  hits_A <- unique(queryHits(findOverlaps(A, B, ignore.strand = TRUE)))
  hits_B <- unique(subjectHits(findOverlaps(A, B, ignore.strand = TRUE)))
  list(
    n_A = length(A),
    n_B = length(B),
    n_A_overlaps_B = length(hits_A),
    n_B_overlaps_A = length(hits_B),
    pct_A_overlaps_B = ifelse(length(A) == 0, NA, 100 * length(hits_A) / length(A)),
    pct_B_overlaps_A = ifelse(length(B) == 0, NA, 100 * length(hits_B) / length(B))
  )
}

overlap_percent_exact <- function(A_ids, B_ids) {
  A_ids <- as.character(A_ids)
  B_ids <- as.character(B_ids)
  inter <- intersect(A_ids, B_ids)
  list(
    n_A = length(A_ids),
    n_B = length(B_ids),
    n_exact_shared = length(inter),
    pct_A_exact = ifelse(length(A_ids) == 0, NA, 100 * length(inter) / length(A_ids)),
    pct_B_exact = ifelse(length(B_ids) == 0, NA, 100 * length(inter) / length(B_ids))
  )
}

# --- Sex vs Sex ---
F_ids <- res_F$region_ids
M_ids <- res_M$region_ids

gr_F <- regionid_to_granges(F_ids, add_chr = TRUE)
gr_M <- regionid_to_granges(M_ids, add_chr = TRUE)

ov_FM_any   <- overlap_percent_any(gr_F, gr_M)
ov_FM_exact <- overlap_percent_exact(F_ids, M_ids)

# --- Load ORIGINAL (combined-sex) ARIMA regions ---
orig_rdata <- "/home/tander10/sternerlab/brain_methylation/Hippocampus_arima_Region_level_SVAdj.RData"
e0 <- new.env()
load(orig_rdata, envir = e0)

if (exists("aging_arima_region", envir = e0)) {
  arima_orig <- get("aging_arima_region", envir = e0)
} else if (exists("aging_arima_region_sx", envir = e0)) {
  arima_orig <- get("aging_arima_region_sx", envir = e0)
} else if (length(ls(envir = e0)) >= 1) {
  arima_orig <- get(ls(envir = e0)[1], envir = e0)
  warning("Using first object in original RData: ", ls(envir = e0)[1])
} else {
  stop("Could not find an ARIMA object in original RData: ", orig_rdata)
}

orig_ids <- rownames(arima_orig$predicted.nonzero)
gr_orig  <- regionid_to_granges(orig_ids, add_chr = TRUE)

# --- Sex vs Original ---
ov_FO_any   <- overlap_percent_any(gr_F, gr_orig)
ov_FO_exact <- overlap_percent_exact(F_ids, orig_ids)

ov_MO_any   <- overlap_percent_any(gr_M, gr_orig)
ov_MO_exact <- overlap_percent_exact(M_ids, orig_ids)

# --- Summarize (EXACT matches) ---
overlap_summary <- bind_rows(
  tibble(
    comparison = "Female vs Male (exact RegionID match)",
    n_A = ov_FM_exact$n_A, n_B = ov_FM_exact$n_B,
    n_shared = ov_FM_exact$n_exact_shared,
    pct_A_shared = ov_FM_exact$pct_A_exact,
    pct_B_shared = ov_FM_exact$pct_B_exact
  ),
  tibble(
    comparison = "Female vs Original (exact RegionID match)",
    n_A = ov_FO_exact$n_A, n_B = ov_FO_exact$n_B,
    n_shared = ov_FO_exact$n_exact_shared,
    pct_A_shared = ov_FO_exact$pct_A_exact,
    pct_B_shared = ov_FO_exact$pct_B_exact
  ),
  tibble(
    comparison = "Male vs Original (exact RegionID match)",
    n_A = ov_MO_exact$n_A, n_B = ov_MO_exact$n_B,
    n_shared = ov_MO_exact$n_exact_shared,
    pct_A_shared = ov_MO_exact$pct_A_exact,
    pct_B_shared = ov_MO_exact$pct_B_exact
  )
)

out_overlap_csv <- file.path(output_root, "Overlap_summary_sex_vs_sex_and_original.csv")
write_csv(overlap_summary, out_overlap_csv)
cat("✅ Wrote overlap summary:", out_overlap_csv, "\n")
print(overlap_summary)

# (Optional) Also save ANY-overlap summary
overlap_any_summary <- bind_rows(
  tibble(
    comparison = "Female vs Male (ANY bp overlap)",
    n_A = ov_FM_any$n_A, n_B = ov_FM_any$n_B,
    n_A_overlaps_B = ov_FM_any$n_A_overlaps_B,
    n_B_overlaps_A = ov_FM_any$n_B_overlaps_A,
    pct_A_overlaps_B = ov_FM_any$pct_A_overlaps_B,
    pct_B_overlaps_A = ov_FM_any$pct_B_overlaps_A
  ),
  tibble(
    comparison = "Female vs Original (ANY bp overlap)",
    n_A = ov_FO_any$n_A, n_B = ov_FO_any$n_B,
    n_A_overlaps_B = ov_FO_any$n_A_overlaps_B,
    n_B_overlaps_A = ov_FO_any$n_B_overlaps_A,
    pct_A_overlaps_B = ov_FO_any$pct_A_overlaps_B,
    pct_B_overlaps_A = ov_FO_any$pct_B_overlaps_A
  ),
  tibble(
    comparison = "Male vs Original (ANY bp overlap)",
    n_A = ov_MO_any$n_A, n_B = ov_MO_any$n_B,
    n_A_overlaps_B = ov_MO_any$n_A_overlaps_B,
    n_B_overlaps_A = ov_MO_any$n_B_overlaps_A,
    pct_A_overlaps_B = ov_MO_any$pct_A_overlaps_B,
    pct_B_overlaps_A = ov_MO_any$pct_B_overlaps_A
  )
)

out_overlap_any_csv <- file.path(output_root, "Overlap_summary_ANYbp_sex_vs_sex_and_original.csv")
write_csv(overlap_any_summary, out_overlap_any_csv)
cat("✅ Wrote ANY-overlap summary:", out_overlap_any_csv, "\n")
print(overlap_any_summary)

# ============================================================
# POST-CLUSTER: BED EXPORT + ChIPseeker ANNOTATION (SEX + CLUSTER)
# ============================================================

library(ChIPseeker)
library(org.Mmu.eg.db)
library(GenomicRanges)
library(GenomeInfoDb)
library(txdbmaker)

# -------------------------------
# Helper functions
# -------------------------------

regionid_to_df <- function(region_ids, add_chr = TRUE) {
  parts <- do.call(rbind, strsplit(as.character(region_ids), "_"))
  df <- as.data.frame(parts, stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end")
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  if (add_chr && !grepl("^chr", df$chr[1])) {
    df$chr <- paste0("chr", df$chr)
  }
  df
}

regionid_to_gr <- function(region_ids, txdb) {
  df <- regionid_to_df(region_ids, add_chr = TRUE)
  gr <- GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end)
  )
  seqlevelsStyle(gr) <- seqlevelsStyle(txdb)
  gr <- keepSeqlevels(
    gr,
    intersect(seqlevels(gr), seqlevels(txdb)),
    pruning.mode = "coarse"
  )
  gr
}

write_region_bed <- function(region_ids, outfile) {
  bed <- regionid_to_df(region_ids, add_chr = TRUE)
  write.table(
    bed,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# -------------------------------
# Build TxDb (once)
# -------------------------------

cat("🔧 Building TxDb (Macaca mulatta)...\n")

txdb <- txdbmaker::makeTxDbFromEnsembl(
  organism = "Macaca mulatta",
  server   = "ensembldb.ensembl.org"
)

# -------------------------------
# Sex-level annotation + BED
# -------------------------------

annotate_and_save_sex <- function(res, sex_label, output_root) {

  outdir <- file.path(output_root, sex_label, "Annotation")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  region_ids <- res$region_ids

  # BED
  bed_file <- file.path(outdir, paste0(sex_label, "_AllRegions.bed"))
  write_region_bed(region_ids, bed_file)

  # Annotation
  gr <- regionid_to_gr(region_ids, txdb)
  annot <- annotatePeak(
    gr,
    TxDb = txdb,
    tssRegion = c(-3000, 3000),
    annoDb = "org.Mmu.eg.db"
  )

  saveRDS(
    annot,
    file.path(outdir, paste0("Annotated_", sex_label, "_AllRegions.rds"))
  )

  cat("✅", sex_label, ":",
      length(region_ids), "regions annotated & BED saved\n")

  invisible(annot)
}

annot_F <- annotate_and_save_sex(res_F, "Female", output_root)
annot_M <- annotate_and_save_sex(res_M, "Male",   output_root)

# -------------------------------
# Cluster-level annotation + BED
# -------------------------------

annotate_clusters <- function(res, sex_label, output_root) {

  clusters <- res$clusters
  region_ids <- names(clusters)

  base_dir <- file.path(output_root, sex_label, "Annotation")

  for (k in sort(unique(clusters))) {

    regs <- region_ids[clusters == k]
    if (length(regs) < 5) next

    cluster_dir <- file.path(base_dir, paste0("Cluster_", k))
    dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)

    # BED
    bed_file <- file.path(
      cluster_dir,
      paste0(sex_label, "_Cluster", k, ".bed")
    )
    write_region_bed(regs, bed_file)

    # Annotation
    gr <- regionid_to_gr(regs, txdb)
    annot <- annotatePeak(
      gr,
      TxDb = txdb,
      tssRegion = c(-3000, 3000),
      annoDb = "org.Mmu.eg.db"
    )

    saveRDS(
      annot,
      file.path(
        cluster_dir,
        paste0("Annotated_", sex_label, "_Cluster", k, ".rds")
      )
    )

    cat("   ↳", sex_label, "Cluster", k, ":",
        length(regs), "regions\n")
  }
}

annotate_clusters(res_F, "Female", output_root)
annotate_clusters(res_M, "Male",   output_root)

cat("✅ All BED files and ChIPseeker annotations completed.\n")
