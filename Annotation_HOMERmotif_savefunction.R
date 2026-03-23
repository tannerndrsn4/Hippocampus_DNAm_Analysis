#Saving annotation files from Chipseeker analysis
# Extract annotated peaks
annotated_peaks <- as.data.frame(VMPs_Decrease@anno)

# Export the main annotated peaks to a CSV file
write.csv(annotated_peaks, 
          file = "VMPs_Decrease_anno.csv", 
          row.names = FALSE)

# Save additional components if they exist in the csAnno object
if (!is.null(VMPs_Decrease@annoStat)) {
  summary_stats <- as.data.frame(VMPs_Decrease@annoStat)
  write.csv(summary_stats, 
            file = "VMPs_Decrease_summary_stats.csv", 
            row.names = FALSE)
}

annotated_peaks <- as.data.frame(VMPs_Increase@anno)

# Export the main annotated peaks to a CSV file
write.csv(annotated_peaks, 
          file = "VMPs_Increase_anno.csv", 
          row.names = FALSE)

# Save additional components if they exist in the csAnno object
if (!is.null(VMPs_Increase@annoStat)) {
  summary_stats <- as.data.frame(VMPs_Increase@annoStat)
  write.csv(summary_stats, 
            file = "VMPs_Increase_summary_stats.csv", 
            row.names = FALSE)
}




#### Function for filtering HOMER motif analysis files to retain only significant and mammalian motifs 
filter_mammalian_motifs <- function(input_file, output_dir = ".", fdr_threshold = 0.01) {
  # Load libraries
  library(dplyr)
  library(stringr)
  
  # Read HOMER knownResults.txt file
  df <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  # Step 1: Remove any motifs corresponding to non-mammalian transcription factors
  non_mammalian_patterns <- c(
    # Plant (Arabidopsis, maize, etc.)
    "Arabidopsis", "col-", "colamp-", "AT[1-9]G", "At[1-9]g", "GSE60143",
    "MYB[0-9]+", "GATA15", "FHY3", "FAR1", "SPL[0-9]+", "RKD", "BBX", "TCP",
    "OBP", "BPC", "CDF", "CAMTA", "COG", "HY5", "O2\\(bZIP\\)", "Corn",
    
    # Yeast / fungi
    "Yeast", "SacCer", "ABF1", "Cbf1",
    
    # Worm / C. elegans
    "cElegans", "modEncode",
    
    # Fly / Drosophila
    "Drosophila", "Embryo-", "WholeFly", "Zelda", "DREF", "SUT1?",
    
    # Amphibian / avian (non-mammalian vertebrates)
    "Frog", "Chicken", "colamp-Frog", "colamp-Chicken"
  )
  
  filtered_df <- df %>%
    filter(
      !str_detect(`Motif Name`, paste(non_mammalian_patterns, collapse = "|")),
      !str_detect(`Motif Name`, "\\(AP2EREBP\\)|\\(NAC\\)|\\(WRKY\\)")  # Keep your earlier filters too
    )
  
  # Step 2: Keep only motifs below FDR threshold
  filtered_df <- filtered_df %>%
    filter(`q-value (Benjamini)` < fdr_threshold)
  
  # Step 3: Count remaining motifs
  n_motifs <- nrow(filtered_df)
  
  # Step 4: Write filtered table to new file
  output_file <- file.path(output_dir, paste0("mammalian_filtered_top", n_motifs, ".txt"))
  write.table(filtered_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message(paste0("✅ Filtered file saved as: ", output_file, " (", n_motifs, " motifs retained)"))
  
  return(filtered_df)
}

filtered_results <- filter_mammalian_motifs("knownResults.txt")