# -------------------------------
# 1. Load Libraries
# -------------------------------
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(AnnotationHub)
library(tidyverse)
library(plyranges)

# -------------------------------
# 2. Load Your Data
# -------------------------------
regions <- read_csv("DNAm_RNAseq_Intersection_distal_intergenic_all_clusters.csv")

gr <- GRanges(seqnames = regions$chr,
              ranges = IRanges(start = regions$start, end = regions$end),
              strand = regions$strand,
              DNAm_Cluster = regions$DNAm_Cluster)

# -------------------------------
# 3. Import UCSC RepeatMasker for Mmul_10
# -------------------------------
# Option A: Download manually from UCSC Table Browser:
#   Group: Repeats
#   Track: RepeatMasker
#   Table: rmsk
#   Assembly: rheMac10 (Mmul_10)
# Save as "rheMac10_rmsk.bed"

repeatmasker <- import("rmsk_mmul10.bed")

# -------------------------------
# 4. Find Overlaps with RepeatMasker
# -------------------------------
hits <- findOverlaps(gr, repeatmasker)

overlap_tbl <- tibble(
  region_chr   = as.character(seqnames(gr)[queryHits(hits)]),
  region_start = as.numeric(start(gr)[queryHits(hits)]),
  region_end   = as.numeric(end(gr)[queryHits(hits)]),
  DNAm_Cluster = as.character(mcols(gr)$DNAm_Cluster[queryHits(hits)]),
  TE_chr       = as.character(seqnames(repeatmasker)[subjectHits(hits)]),
  TE_start     = as.numeric(start(repeatmasker)[subjectHits(hits)]),
  TE_end       = as.numeric(end(repeatmasker)[subjectHits(hits)]),
  TE_name      = as.character(mcols(repeatmasker)$name[subjectHits(hits)])
)

# Collapse multiple TEs per region
te_overlap_collapsed <- overlap_tbl %>%
  group_by(region_chr, region_start, region_end, DNAm_Cluster) %>%
  summarise(
    TE_chrs     = paste(unique(na.omit(TE_chr)), collapse = "; "),
    TE_starts   = paste(unique(na.omit(TE_start)), collapse = "; "),
    TE_ends     = paste(unique(na.omit(TE_end)), collapse = "; "),
    TE_names    = paste(unique(na.omit(TE_name)), collapse = "; "),
  ) %>%
  ungroup()

dim(overlap_tbl)

write_csv(overlap_tbl, "Distal_intergenic_TE_overlaps.csv")

# -------------------------------
# 5. Enrichment Test for TE Overlap
# -------------------------------

library(regioneR)

# Your objects:
# gr            = your distal intergenic query regions (GRanges)
# repeatmasker  = TE annotation (GRanges)
# background    = your annotated_background.rds object

background <- readRDS("annotated_background.rds")
# Extract the GRanges from background
background_gr <- background@anno  # this is the GRanges we’ll use as the universe

# Define genome
genome_mmul10 <- getGenomeAndMask("rheMac10")

# Run permutation test with custom background
set.seed(123)
perm_test_TE <- overlapPermTest(
  A = gr,
  B = repeatmasker,
  genome = genome_mmul10$genome,
  universe = background_gr,   # restrict randomization to your background regions
  ntimes = 100,
  alternative = "greater"
)

# View result
perm_test_TE
pval <- perm_test_TE$numOverlaps$pval
pval

plot(perm_test_TE)

# -------------------------------
# 6. LiftOver to Human (Mmul_10 → hg38)
# -------------------------------
chain <- import.chain("rheMac10ToHg38.over.chain")

# --- Ensure chromosome names have "chr" prefix ---
# Fix for both query (gr) and background (background@anno)
seqlevelsStyle(gr) <- "UCSC"
seqlevelsStyle(background@anno) <- "UCSC"

# --- Lift over the distal intergenic query regions ---
gr_hg38_list <- liftOver(gr, chain)
gr_hg38 <- unlist(gr_hg38_list[elementNROWS(gr_hg38_list) == 1])

# --- Lift over the background regions ---
background_list <- liftOver(background@anno, chain)
background_hg38 <- unlist(background_list[elementNROWS(background_list) == 1])

cat("✅ LiftOver complete:",
    length(gr_hg38), "query regions and",
    length(background_hg38), "background regions successfully mapped to hg38\n")

# -------------------------------
# 7. Import ENCODE Enhancer Regions (hg38)
# -------------------------------
# Read manually (it's tab-separated)
enh_df <- read_tsv("GRCh38-ELS.bed", 
                   col_names = c("chr", "start", "end", "dataset_ID", "enhancer_ID", "enhancer_type"),
                   col_types = cols(
                     chr = col_character(),
                     start = col_double(),
                     end = col_double(),
                     dataset_ID = col_character(),
                     enhancer_ID = col_character(),
                     enhancer_type = col_character()
                   ))

# Convert to GRanges
encode_enh <- makeGRangesFromDataFrame(
  enh_df,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# -------------------------------
# 8. Find Overlaps with ENCODE Enhancers
# -------------------------------
enh_hits <- findOverlaps(gr_hg38, encode_enh)

enh_overlap_tbl <- tibble(
  region_chr   = as.character(seqnames(gr_hg38)[queryHits(enh_hits)]),
  region_start = as.numeric(start(gr_hg38)[queryHits(enh_hits)]),
  region_end   = as.numeric(end(gr_hg38)[queryHits(enh_hits)]),
  DNAm_Cluster = as.character(mcols(gr_hg38)$DNAm_Cluster[queryHits(enh_hits)]),  # <-- changed here
  
  Enhancer_chr   = as.character(seqnames(encode_enh)[subjectHits(enh_hits)]),
  Enhancer_start = as.numeric(start(encode_enh)[subjectHits(enh_hits)]),
  Enhancer_end   = as.numeric(end(encode_enh)[subjectHits(enh_hits)]),
  Enhancer_dataset_ID = as.character(mcols(encode_enh)$dataset_ID[subjectHits(enh_hits)]),
  Enhancer_ID         = as.character(mcols(encode_enh)$enhancer_ID[subjectHits(enh_hits)]),
  Enhancer_type       = as.character(mcols(encode_enh)$enhancer_type[subjectHits(enh_hits)])
)

# --- 3D chromatin links ---
chromatin_links <- read_tsv("V4-hg38.Gene-Links.3D-Chromatin.txt",
                            col_names = c("Enhancer_ID", "gene_id", "gene_name", "gene_type",
                                          "method", "experiment_ID", "cell_type", "link_score", "p_value"),
                            col_types = cols(.default = col_character())
) %>%
  mutate(link_source = "3D_Chromatin")

# --- CRISPR links ---
crispr_links <- read_tsv("V4-hg38.Gene-Links.CRISPR.txt",
                         col_names = c("Enhancer_ID", "gene_id", "gene_name", "gene_type",
                                       "coord", "method", "experiment_ID", "cell_type", "score", "p_value"),
                         col_types = cols(.default = col_character())
) %>%
  mutate(link_source = "CRISPR")

# --- eQTL links ---
eqtl_links <- read_tsv("V4-hg38.Gene-Links.eQTLs.txt",
                       col_names = c("Enhancer_ID", "gene_id", "gene_name", "gene_type",
                                     "variant", "method", "tissue", "effect_size", "p_value"),
                       col_types = cols(.default = col_character())
) %>%
  mutate(link_source = "eQTL")

# Combine them all into one big table
all_links <- bind_rows(chromatin_links, crispr_links, eqtl_links)

enh_overlap_annot <- enh_overlap_tbl %>%
  left_join(all_links, by = "Enhancer_ID")

enh_overlap_collapsed <- enh_overlap_annot %>%
  group_by(region_chr, region_start, region_end, DNAm_Cluster,
           Enhancer_chr, Enhancer_start, Enhancer_end,
           Enhancer_dataset_ID, Enhancer_ID, Enhancer_type) %>%
  summarise(
    linked_genes = paste(unique(na.omit(gene_name)), collapse = "; "),
    link_sources = paste(unique(na.omit(link_source)), collapse = "; ")
  ) %>%
  ungroup()

# Save results
write_csv(enh_overlap_collapsed, "LiftedOver_ENCODE_Enhancer_Overlaps.csv")

# Count unique regions in TE overlap table
num_unique_TE_regions <- overlap_tbl %>%
  distinct(region_chr, region_start, region_end) %>%
  nrow()

# Count unique regions in enhancer overlap table
num_unique_enh_regions <- enh_overlap_collapsed %>%
  distinct(region_chr, region_start, region_end) %>%
  nrow()

# Print the results
cat("Unique regions in TE overlap table:", num_unique_TE_regions, "\n")
cat("Unique regions in enhancer overlap table:", num_unique_enh_regions, "\n")

# -------------------------------
# 9. Enrichment Test for Enhancer Overlap (using lifted-over custom background)
# -------------------------------
set.seed(123)

# Restrict all to common sequence levels
common_chrs <- intersect(seqlevels(gr_hg38),
                         intersect(seqlevels(background_hg38), seqlevels(encode_enh)))

gr_hg38 <- keepSeqlevels(gr_hg38, common_chrs, pruning.mode = "coarse")
background_hg38 <- keepSeqlevels(background_hg38, common_chrs, pruning.mode = "coarse")
encode_enh <- keepSeqlevels(encode_enh, common_chrs, pruning.mode = "coarse")

# Define genome (let regioneR infer lengths automatically)
genome_from_bg <- getGenomeAndMask("hg38")$genome

# Run permutation test
perm_test_ENC <- overlapPermTest(
  A = gr_hg38,                  # your lifted-over query regions
  B = encode_enh,               # ENCODE enhancers
  universe = background_hg38,   # custom background lifted over to hg38
  genome = genome_from_bg,
  ntimes = 1000,
  alternative = "greater"
)

# Extract and display results
pval_enh <- perm_test_ENC$numOverlaps$pval
zscore_enh <- perm_test_ENC$numOverlaps$zscore

cat("Observed overlaps:", perm_test_ENC$numOverlaps$observed, "\n")
cat("Mean random overlaps:", mean(perm_test_ENC$numOverlaps$permuted), "\n")
cat("Z-score:", zscore_enh, "\n")
cat("Permutation p-value:", pval_enh, "\n")

if (pval_enh < 0.05) {
  message("✅ Significant enrichment for ENCODE enhancers (p < 0.05)")
} else {
  message("❌ No significant enrichment for ENCODE enhancers.")
}

# -------------------------------
# 10. Identify Regions Overlapping Both TEs and Enhancers
# -------------------------------

cat("\n🔍 Checking for regions overlapping both TEs and enhancers...\n")

# --- Step 10.1: Map lifted-over enhancer regions back to macaque coordinates ---
# (We invert the chain to go hg38 → rheMac10)
chain_rev <- import.chain("hg38ToRheMac10.over.chain")

# Lift enhancer-overlapping regions back
enh_liftback_list <- liftOver(
  GRanges(seqnames = enh_overlap_collapsed$region_chr,
          ranges   = IRanges(start = enh_overlap_collapsed$region_start,
                             end   = enh_overlap_collapsed$region_end)),
  chain_rev
)

# Keep only one-to-one mapped regions for clean comparison
enh_liftback <- unlist(enh_liftback_list[elementNROWS(enh_liftback_list) == 1])
cat("Lifted back", length(enh_liftback), "enhancer-overlapping regions to macaque coordinates.\n")

# --- Step 10.2: Convert TE-overlapping regions to GRanges ---
te_gr <- GRanges(seqnames = overlap_tbl$region_chr,
                 ranges   = IRanges(start = overlap_tbl$region_start,
                                    end   = overlap_tbl$region_end))

# --- Step 10.3: Find overlaps between TE and enhancer regions (in macaque coords) ---
both_hits <- findOverlaps(te_gr, enh_liftback)

num_shared <- length(unique(queryHits(both_hits)))
num_TE <- length(unique(te_gr))
num_enh <- length(unique(enh_liftback))

cat("✅ Regions overlapping both TEs and enhancers:", num_shared, "\n")
cat("   → Fraction of TE regions also enhancer-linked:",
    round(num_shared / num_TE, 4), "\n")
cat("   → Fraction of enhancer-linked regions also TE-linked:",
    round(num_shared / num_enh, 4), "\n")

# --- Optional: Save list of shared regions ---
shared_regions <- tibble(
  chr   = as.character(seqnames(te_gr)[queryHits(both_hits)]),
  start = start(te_gr)[queryHits(both_hits)],
  end   = end(te_gr)[queryHits(both_hits)]
)

write_csv(shared_regions, "Shared_TE_and_Enhancer_Regions_rheMac10.csv")

cat("💾 Saved shared region list: Shared_TE_and_Enhancer_Regions_rheMac10.csv\n")