#Gprofiler on bed files 
# Load required packages
library(gprofiler2)
library(GenomicRanges)
library(biomaRt)
library(dplyr)
library(stringr)

# ---- STEP 1: Load BED file ----
bed <- read.delim("AllVMPs_M_HOMER.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(bed) <- c("chr", "start", "end", "name", "score", "strand")

# ---- STEP 2: Fix chromosome names (remove 'chr' prefix) ----
bed$chr <- str_replace(bed$chr, "^chr", "")

# ---- STEP 3: Create GRanges object ----
gr <- GRanges(seqnames = bed$chr,
              ranges = IRanges(start = bed$start, end = bed$end))

# ---- STEP 4: Map CpGs to nearest genes using biomaRt (Macaca mulatta genome) ----
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmulatta_gene_ensembl")

annot <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "external_gene_name"),
  mart = mart
)

# Convert to GRanges
annot_gr <- GRanges(seqnames = annot$chromosome_name,
                    ranges = IRanges(start = annot$start_position, end = annot$end_position),
                    gene_id = annot$ensembl_gene_id,
                    gene_name = annot$external_gene_name)

# ---- STEP 5: Align seqlevels between gr and annot_gr ----
common_seqlevels <- intersect(seqlevels(gr), seqlevels(annot_gr))
gr <- keepSeqlevels(gr, common_seqlevels, pruning.mode = "coarse")
annot_gr <- keepSeqlevels(annot_gr, common_seqlevels, pruning.mode = "coarse")

# ---- STEP 6: Find overlaps between your BED regions and genes ----
hits <- findOverlaps(gr, annot_gr)
genes <- unique(mcols(annot_gr)$gene_name[subjectHits(hits)])

cat("✅ Found", length(genes), "unique overlapping genes.\n")

# ---- STEP 7: Run g:Profiler enrichment ----
if (length(genes) > 0) {
  gostres <- gost(
    query = genes,
    organism = "mmulatta",
    sources = c("GO:BP", "KEGG"),
    correction_method = "fdr",
    significant = TRUE,
    user_threshold = 0.05,
    domain_scope = "annotated"
  )
  
  # ---- STEP 8: Filter and unlist before saving ----
  if (!is.null(gostres$result)) {
    gostres_df <- gostres$result %>%
      filter(term_size <= 2500)
    
    # Convert list columns to comma-separated strings
    gostres_df[] <- lapply(gostres_df, function(x) {
      if (is.list(x)) sapply(x, paste, collapse = ",") else x
    })
    
    # Save results as CSV
    write.csv(gostres_df, "gProfiler_AllVMPs_M.csv", row.names = FALSE)
    
    cat("✅ Enrichment analysis complete —", nrow(gostres_df),
        "terms passed FDR < 0.05 and term_size ≤ 2500.\nResults saved to gProfiler_VMPs_increase.csv\n")
  } else {
    cat("⚠️ No significant enrichment terms found (FDR < 0.05).\n")
  }
} else {
  cat("⚠️ No overlapping genes found between BED regions and Macaca mulatta annotation.\n")
}