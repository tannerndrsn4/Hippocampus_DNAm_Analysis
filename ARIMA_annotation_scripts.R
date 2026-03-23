# ============================================================
# Script: VMP_annotation_summary.R
# Purpose: Combine ARIMA cluster annotations and create summary + BED files
# ============================================================

# --- Load libraries ---
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# --- User Inputs ---
input_rdata <- "Annotated_REGION_Clusters.RData"
output_prefix <- "master_VMP"

# --- Load annotated cluster file ---
cat("📂 Loading annotation data from", input_rdata, "...\n")
load(input_rdata)

# --- Extract the @anno data frame from each cluster ---
# Each annotated_clusters$ClusterX is a csAnno object
cat("✅ Extracting annotation tables from 4 clusters...\n")
cluster1 <- as.data.frame(annotated_clusters$Cluster1@anno)
cluster2 <- as.data.frame(annotated_clusters$Cluster2@anno)
cluster3 <- as.data.frame(annotated_clusters$Cluster3@anno)
cluster4 <- as.data.frame(annotated_clusters$Cluster4@anno)

# --- Flatten list columns to character ---
flatten_df <- function(df) {
  df[] <- map(df, function(x) {
    if (is.list(x)) {
      sapply(x, function(i) paste(i, collapse = "; "))
    } else {
      x
    }
  })
  df
}

cluster1 <- flatten_df(cluster1)
cluster2 <- flatten_df(cluster2)
cluster3 <- flatten_df(cluster3)
cluster4 <- flatten_df(cluster4)

# --- Combine all clusters ---
master_VMP <- bind_rows(cluster1, cluster2, cluster3, cluster4)
cat("✅ Combined all clusters (total rows:", nrow(master_VMP), ")\n")

# --- Ensure chromosome names are prefixed with 'chr' ---
master_VMP <- master_VMP %>%
  mutate(chr = ifelse(grepl("^chr", seqnames), seqnames, paste0("chr", seqnames)))

# --- Extract simplified feature types from annotation ---
extract_feature <- function(anno) {
  case_when(
    str_detect(anno, regex("Promoter", ignore_case = TRUE)) ~ "Promoter",
    str_detect(anno, regex("Intron", ignore_case = TRUE)) ~ "Intron",
    str_detect(anno, regex("Exon", ignore_case = TRUE)) ~ "Exon",
    str_detect(anno, regex("UTR", ignore_case = TRUE)) ~ "UTR",
    str_detect(anno, regex("Downstream", ignore_case = TRUE)) ~ "Downstream",
    str_detect(anno, regex("Intergenic", ignore_case = TRUE)) ~ "Distal Intergenic",
    TRUE ~ "Other"
  )
}

master_VMP <- master_VMP %>%
  mutate(feature = extract_feature(annotation))

# --- Save master annotation file ---
write_csv(master_VMP, paste0(output_prefix, "_annotation.csv"))
cat("✅ Combined master annotation file saved as", paste0(output_prefix, "_annotation.csv"), "\n")

# --- Summarize feature frequencies ---
annotation_summary <- master_VMP %>%
  mutate(feature = as.character(feature)) %>%
  group_by(feature) %>%
  summarise(Count = n()) %>%
  mutate(Frequency = 100 * Count / sum(Count)) %>%
  arrange(desc(Frequency))

# --- Save summary ---
write_csv(annotation_summary, paste0(output_prefix, "_annotation_summary.csv"))
cat("✅ Genomic annotation summary saved as", paste0(output_prefix, "_annotation_summary.csv"), "\n\n")

print(annotation_summary)

# --- Generate feature-specific BED files ---
for (f in unique(master_VMP$feature)) {
  subset_df <- master_VMP %>%
    filter(feature == f) %>%
    select(chr, start, end)
  
  if (nrow(subset_df) > 0) {
    bed_file <- paste0(tolower(gsub(" ", "_", f)), "_regions.bed")
    write_tsv(subset_df, bed_file, col_names = FALSE)
    cat("✅ BED file created:", bed_file, "(", nrow(subset_df), "entries)\n")
  }
}

cat("\n🎉 All done! BED files and summary stats have been generated successfully.\n")