library(dplyr)
library(readr)
library(stringr)
library(purrr)

input_rdata <- "Annotated_REGION_Clusters.RData"
deg_file <- "RNAseq_genes.csv"
output_prefix <- "DNAm_RNAseq_Intersection"

default_rna_cluster <- "RNAseq_Cluster1"

# ============================================================
# Step 1: Load annotated clusters
# ============================================================
cat("📂 Loading annotated clusters from", input_rdata, "...\n")
load(input_rdata)

# --- Helper: flatten list-columns into simple character values ---
flatten_df <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.list(x)) {
      sapply(x, function(i) paste(i, collapse = "; "), USE.NAMES = FALSE)
    } else {
      as.character(x)
    }
  })
  as.data.frame(df)
}

# --- Extract and flatten annotation tables from each cluster ---
cat("✅ Extracting and flattening annotation tables for 4 clusters...\n")
cluster1 <- flatten_df(as.data.frame(annotated_clusters$Cluster1@anno))
cluster2 <- flatten_df(as.data.frame(annotated_clusters$Cluster2@anno))
cluster3 <- flatten_df(as.data.frame(annotated_clusters$Cluster3@anno))
cluster4 <- flatten_df(as.data.frame(annotated_clusters$Cluster4@anno))

# --- Add DNAm cluster labels before merging ---
cluster1$DNAm_Cluster <- "Cluster1"
cluster2$DNAm_Cluster <- "Cluster2"
cluster3$DNAm_Cluster <- "Cluster3"
cluster4$DNAm_Cluster <- "Cluster4"

# --- Combine all clusters ---
master_annot <- bind_rows(cluster1, cluster2, cluster3, cluster4)
cat("✅ Combined all clusters into master annotation (", nrow(master_annot), " regions )\n")

# --- Ensure chromosome names have 'chr' prefix ---
master_annot <- master_annot %>%
  mutate(chr = ifelse(grepl("^chr", seqnames), seqnames, paste0("chr", seqnames)))

# ============================================================
# Step 2: Load RNAseq DEGs
# ============================================================
deg_list <- read_csv(deg_file, show_col_types = FALSE)

# Determine RNAseq cluster info
if (ncol(deg_list) >= 2) {
  # If RNAseq cluster info is provided
  names(deg_list)[1:2] <- c("SYMBOL", "RNAseq_Cluster")
  deg_symbols <- unique(deg_list$SYMBOL)
  cat("✅ Loaded", length(deg_symbols), "unique DEGs across",
      length(unique(deg_list$RNAseq_Cluster)), "RNAseq clusters from", deg_file, "\n")
} else {
  # If only gene symbols provided
  deg_list <- deg_list %>%
    rename(SYMBOL = 1) %>%
    mutate(RNAseq_Cluster = default_rna_cluster)
  deg_symbols <- unique(deg_list$SYMBOL)
  cat("✅ Loaded", length(deg_symbols), "unique DEGs (assigned to default cluster:",
      default_rna_cluster, ")\n")
}

# ============================================================
# Step 3: Intersect DNAm annotations with DEG symbols
# ============================================================
master_annot$SYMBOL <- as.character(master_annot$SYMBOL)

candidate_regions <- master_annot %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(SYMBOL %in% deg_symbols)

cat("✅ Found", nrow(candidate_regions),
    "candidate DNAm regions overlapping DEGs (after removing NA SYMBOLs)\n")

# Merge RNAseq cluster info onto matching genes
candidate_regions <- candidate_regions %>%
  left_join(deg_list %>% select(SYMBOL, RNAseq_Cluster), by = "SYMBOL")

# ============================================================
# Step 4: Count duplicates (DEGs with multiple regions)
# ============================================================
duplicate_genes <- candidate_regions %>%
  group_by(SYMBOL) %>%
  summarise(NumRegions = n()) %>%
  filter(NumRegions > 1)

cat("🔁 There are", nrow(duplicate_genes),
    "DEGs with multiple associated DNAm regions\n")

# ============================================================
# Step 5: Extract simplified feature categories
# ============================================================
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

candidate_regions <- candidate_regions %>%
  mutate(feature = extract_feature(annotation))

# ============================================================
# Step 6: Summarize by feature
# ============================================================
feature_summary <- candidate_regions %>%
  group_by(feature) %>%
  summarise(Count = n()) %>%
  mutate(Frequency = 100 * Count / sum(Count)) %>%
  arrange(desc(Frequency))

cat("\n📊 Candidate regulatory region feature breakdown:\n")
print(feature_summary)

# ============================================================
# Step 6b: Subset to genic features only (remove intergenic)
# ============================================================

genic_features <- c("Promoter", "Intron", "Exon", "UTR")

candidate_regions_genic <- candidate_regions %>%
  filter(feature %in% genic_features)

cat("🧬 Retained", nrow(candidate_regions_genic),
    "genic candidate regions (excluded distal intergenic and downstream)\n")

# ============================================================
# Step 6c: Count duplicate genes (genic-only)
# ============================================================

duplicate_genes_genic <- candidate_regions_genic %>%
  group_by(SYMBOL) %>%
  summarise(NumRegions = n()) %>%
  filter(NumRegions > 1)

cat("🔁 There are", nrow(duplicate_genes_genic),
    "DEGs with multiple genic DNAm regions\n")

# ============================================================
# Step 7: Save outputs
# ============================================================
write_csv(candidate_regions, paste0(output_prefix, "_candidates.csv"))
write_csv(feature_summary, paste0(output_prefix, "_feature_summary.csv"))

cat("\n💾 Candidate regulatory regions saved as:",
    paste0(output_prefix, "_candidates.csv"), "\n")
cat("💾 Feature summary saved as:",
    paste0(output_prefix, "_feature_summary.csv"), "\n")

write_csv(candidate_regions_genic,
          paste0(output_prefix, "_candidates_genic_only.csv"))

write_csv(duplicate_genes_genic,
          paste0(output_prefix, "_duplicate_genes_genic_only.csv"))

cat("💾 Genic-only candidate regions saved as:",
    paste0(output_prefix, "_candidates_genic_only.csv"), "\n")

cat("💾 Genic-only duplicate gene summary saved as:",
    paste0(output_prefix, "_duplicate_genes_genic_only.csv"), "\n")

# Reuse the same feature extraction logic you already defined later
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

distal_intergenic_regions <- master_annot %>%
  mutate(feature = extract_feature(annotation)) %>%
  filter(feature == "Distal Intergenic") %>%
  mutate(
    Region_Type = "distal intergenic"  # column that is the same for all rows
  ) %>%
  # keep what you asked for: region location, cluster, and the label column
  select(chr, start, end, DNAm_Cluster, Region_Type, everything())

distal_outfile <- paste0(output_prefix, "_distal_intergenic_all_clusters.csv")
write_csv(distal_intergenic_regions, distal_outfile)

cat("💾 Distal intergenic regions saved as:", distal_outfile,
    "(", nrow(distal_intergenic_regions), "regions )\n")

# ============================================================
# Step 8: Save feature-specific BED files
# ============================================================
for (f in unique(candidate_regions$feature)) {
  subset_df <- candidate_regions %>%
    filter(feature == f) %>%
    select(chr, start, end)
  
  if (nrow(subset_df) > 0) {
    bed_file <- paste0(tolower(gsub(" ", "_", f)), "_candidate_regions.bed")
    write_tsv(subset_df, bed_file, col_names = FALSE)
    cat("✅ BED file created:", bed_file, "(", nrow(subset_df), "entries)\n")
  }
}

cat("\n🎉 Intersection analysis complete! Candidate region tables and feature-specific BEDs generated.\n")
