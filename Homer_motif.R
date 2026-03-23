# Set working directory to where your BED files are
setwd("/Users/tanneranderson/Desktop/PhD brain evolution/DNA methylation/ARIMA results")

format_bed_for_homer <- function(bedfile, output_file = NULL) {
  bed <- read.table(bedfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end")
  
  # Create unique IDs
  bed$peak_id <- paste(bed$chr, bed$start, bed$end, sep = "_")
  
  # Add dummy score & strand columns
  bed$score <- "."
  bed$strand <- "+"
  
  # Re-order to HOMER’s 6-column format
  bed_out <- bed[, c("chr", "start", "end", "peak_id", "score", "strand")]
  
  if (is.null(output_file))
    output_file <- sub(".bed$", "_homer.bed", bedfile)
  
  write.table(
    bed_out,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  message("Formatted BED saved: ", output_file)
}

# Apply to all your cluster files
files <- list.files(pattern = "Cluster[0-9]_regions\\.bed$")
lapply(files, format_bed_for_homer)

# Optionally also format the background file
format_bed_for_homer("Background_CpGs.bed", "Background_CpGs_homer.bed")

### RUN NEXT PART IN TERMINAL ###
for FILE in *HOMER.bed; do
BASENAME=$(basename "$FILE" _regions_homer.bed)
echo "Running HOMER motif enrichment for $BASENAME (75 bp fixed window)..."

findMotifsGenome.pl \
"$FILE" rheMac10 \
"homer_motif_results/${BASENAME}_motifs" \
-bg "Background_CpGs_homer.bed" \
-size 75 \
-len 8,10,12 \
-p 8

echo "$BASENAME analysis complete."
echo "Results saved to homer_motif_results/${BASENAME}_motifs"
echo "------------------------------------------------------------"
done