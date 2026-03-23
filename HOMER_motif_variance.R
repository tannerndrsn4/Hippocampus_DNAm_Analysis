# Set working directory to the Variance Analysis folder
setwd("/Users/tanneranderson/Desktop/PhD brain evolution/DNA methylation/Variance Analysis")

# Function to reformat .bed files for HOMER
format_bed_for_homer <- function(bedfile, output_file = NULL) {
  bed <- read.table(bedfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(bed) <- c("chr", "start", "end")
  
  # Create unique peak IDs
  bed$peak_id <- paste(bed$chr, bed$start, bed$end, sep = "_")
  
  # Add dummy score & strand columns
  bed$score <- "."
  bed$strand <- "+"
  
  # Reorder to HOMER’s 6-column format
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

# Reformat the variance-related BED files
format_bed_for_homer("VMPs_DecreasingVariance_Age.bed", "VMPs_DecreasingVariance_Age_homer.bed")
format_bed_for_homer("VMPs_IncreasingVariance_Age.bed", "VMPs_IncreasingVariance_Age_homer.bed")

# Reformat your background CpGs BED file
format_bed_for_homer("background_cpgs.bed", "background_cpgs_homer.bed")

cd "/Users/tanneranderson/Desktop/PhD brain evolution/DNA methylation/Variance Analysis"

# Create an output directory for motif results
mkdir -p homer_motif_results

# Loop through the variance BED files and run HOMER
for FILE in VMPs_*Variance_Age_homer.bed; do
BASENAME=$(basename "$FILE" _homer.bed)
echo "Running HOMER motif enrichment for $BASENAME (75 bp fixed window)..."

findMotifsGenome.pl \
"$FILE" rheMac10 \
"homer_motif_results/${BASENAME}_motifs" \
-bg "background_cpgs_homer.bed" \
-size 75 \
-len 8,10,12 \
-p 8

echo "$BASENAME analysis complete."
echo "Results saved to homer_motif_results/${BASENAME}_motifs"
echo "------------------------------------------------------------"
done