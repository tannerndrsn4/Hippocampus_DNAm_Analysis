# Load necessary libraries
library(glmnet)
library(bsseq)
library(GenomicRanges)
library(ggplot2)

# Step 1: Load and Normalize the Training Data (perc_filtered)
perc_filtered <- readRDS(file="/home/tander10/sternerlab/brain_methylation/RRBSperc_filtered.rds")
perc_filtered[is.na(perc_filtered)] <- 0

# Step 2: Load Metadata
liverMeta <- read.csv(file="/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata.csv")

# Step 3: Construct the bsseq object with test data
TestMeta <- read.csv(file="/home/tander10/sternerlab/brain_methylation/Rapa_meta.csv")
mmul10_fa <- readRDS(file="/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")

# Filter for the first 21 chromosomes (chr1-chrX)
names(mmul10_fa)[1:21] -> chr
mmul10_fa <- mmul10_fa[chr]

# Find CG sites in the genome
loci <- findLoci("CG", subject=mmul10_fa)

# Read in test Bismark files
bismarkBSseq_test <- read.bismark(
  files = list.files("/home/tander10/sternerlab/Rapamycin/methyl_extract", full.names = TRUE, pattern=".cov.gz"),
  strandCollapse = FALSE,
  verbose = TRUE,
  BACKEND = "HDF5Array",
  loci = loci,
  rmZeroCov = FALSE,
  dir="/projects/sternerlab/shared/Coverage_files/hdf5",
  replace = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)

saveRDS(bismarkBSseq_test, file="/home/tander10/sternerlab/brain_methylation/bismarkBSseq_test_rapa.rds")

# Step 4: Extract coverage information for the windows of interest in the test data
coordinates <- strsplit(rownames(perc_filtered), "_")
chromosomes <- sapply(coordinates, function(coord) coord[1])
starts <- sapply(coordinates, function(coord) as.numeric(coord[2]))
ends <- sapply(coordinates, function(coord) as.numeric(coord[3]))

# Create the GRanges object for windowed regions
reg <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))

coverage_test <- getCoverage(bismarkBSseq_test, regions = reg, type = "Cov", what = "perRegionTotal")
M_vals_test <- getCoverage(bismarkBSseq_test, regions = reg, type = "M", what = "perRegionTotal")
perc_complete_test <- M_vals_test / coverage_test

# Function to set GRanges as row names of a matrix with modified formatting
set_granges_as_row_names <- function(matrix_data, granges_object) {
  # Extract seqnames and ranges from the GRanges object
  seqnames <- seqnames(granges_object)
  ranges <- ranges(granges_object)
  
  # Create a character vector with modified row names
  row_names <- paste(seqnames, ranges@start, ranges@start+999, sep = "_")
  
  # Set row names of the matrix
  rownames(matrix_data) <- row_names
  
  return(matrix_data)
}

coverage_test <- set_granges_as_row_names(coverage_test, reg)
M_vals_test <- set_granges_as_row_names(M_vals_test, reg)
perc_complete_test <- set_granges_as_row_names(perc_complete_test, reg)
perc_complete_test[is.na(perc_complete_test)] <- 0

saveRDS(perc_complete_test, file="/home/tander10/sternerlab/brain_methylation/perc_complete_test_rapa.rds")

# Step 5: Remove Rows with NAs in Both Training and Test Datasets
# Identify rows with NAs in either dataset
na_rows_train <- apply(perc_filtered, 1, function(x) any(is.na(x)))
na_rows_test <- apply(perc_complete_test, 1, function(x) any(is.na(x)))

# Combine rows with NAs from both datasets
rows_to_remove <- na_rows_train | na_rows_test
num_rows_removed <- sum(rows_to_remove)
cat("Number of rows removed due to NAs:", num_rows_removed, "\n")

# Filter out rows with NAs from both datasets
perc_filtered <- perc_filtered[!rows_to_remove, ]
perc_complete_test <- perc_complete_test[!rows_to_remove, ]

# Step 6: Normalize the Training and Test Data
# Normalize the training data
norm_train <- apply(perc_filtered, 1, function(x) {return(qqnorm(x, plot = FALSE)$x)})
norm_train <- t(norm_train)
norm_train <- as.matrix(norm_train)

norm_train2 <- apply(norm_train, 2, function(x) {return(qqnorm(x, plot = FALSE)$x)})
norm_train2 <- t(norm_train2)
colnames(norm_train2) <- colnames(norm_train)

norm_test<-apply(perc_complete_test,1,function(x){return(qqnorm(x,plot=F)$x)})
norm_test <- t(norm_test)
norm_test <- as.matrix(norm_test)

norm_test2<-as.matrix(apply(norm_test,2,function(x){return(qqnorm(x,plot=F)$x)}))
norm_test2<-t(norm_test2)
colnames(norm_test2) <- colnames(norm_test)

# # Step 7: Align Features in Test and Training Data
# matching_features <- intersect(rownames(norm_train2), rownames(norm_test2))
# norm_train2 <- norm_train2[matching_features, , drop = FALSE]
# norm_test2 <- norm_test2[matching_features, , drop = FALSE]
# 
# # Ensure alignment of feature order
# norm_test2 <- norm_test2[rownames(norm_train2), , drop = FALSE]

# Step 8: Train the Model
finalModel <- cv.glmnet(data.matrix(t(perc_filtered)), y = liverMeta$Age, family = "gaussian", alpha = 0.1)

# Step 9: Predict Ages Using the Model
predicted_ages <- predict(finalModel, t(perc_complete_test), s = "lambda.min")
predicted_ages
saveRDS(predicted_ages, file="/home/tander10/sternerlab/brain_methylation/predicted_ages_rapa.rds")


# Step 10: Create and Save a Data Frame with Predicted Ages
predicted_ages_df <- data.frame(
  Sample_ID = TestMeta$Animal.ID,
  Chronological_Age = TestMeta$Age.End.of.Study,
  Predicted_Age = predicted_ages,
  Group = TestMeta$Group,
  Year = TestMeta$Year
)

saveRDS(predicted_ages_df, file="/home/tander10/sternerlab/brain_methylation/predicted_ages_summary_rapa.rds")

# Step 11: Train the Model
finalModel <- cv.glmnet(data.matrix((norm_train2)), y = liverMeta$Age, family = "gaussian", alpha = 0.1)

# Step 12: Predict Ages Using the Model
predicted_ages <- predict(finalModel, (norm_test2), s = "lambda.min")
predicted_ages
saveRDS(predicted_ages, file="/home/tander10/sternerlab/brain_methylation/predicted_ages_rapa_normalized.rds")


# Step 13: Create and Save a Data Frame with Predicted Ages
predicted_ages_df <- data.frame(
  Sample_ID = TestMeta$Animal.ID,
  Chronological_Age = TestMeta$Age.End.of.Study,
  Predicted_Age = predicted_ages,
  Group = TestMeta$Group,
  Year = TestMeta$Year
)

saveRDS(predicted_ages_df, file="/home/tander10/sternerlab/brain_methylation/predicted_ages_summary_rapa_normalized.rds")

# Step 11: Optional Plot
ggplot(predicted_ages_df, aes(x = Chronological_Age, y = Predicted_Age)) +
  geom_point() +
  stat_smooth(method = "lm", col = "dark green") +
  theme_light()