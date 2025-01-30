#1. Script for loading in your data and creating a bsseq object from which you can generate an epigenetic clock. You will need bismark aligned coverage files "cov.gz" for this step.

library(Biostrings)
library(bsseq)
library(ggplot2)


mmul10_fa <- readRDS(file="/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")

## filter for the first 20 chromosomes (chr1-chrX)
names(mmul10_fa)[1:21]->chr
mmul10_fa=mmul10_fa[chr]

## find CG sites in the genome
loci = findLoci("CG", subject=mmul10_fa)

bismarkBSseq <- read.bismark(files = list.files("/home/tander10/sternerlab/brain_methylation/Bismark_output",
                                                full.names = T, pattern=".cov.gz"),
                             strandCollapse = FALSE,
                             verbose = TRUE,
                             BACKEND = "HDF5Array",
                             loci = loci,
                             rmZeroCov = FALSE,
                             dir="/projects/sternerlab/shared/Coverage_files/hdf5",
                             replace = TRUE,
                             BPPARAM = BiocParallel::SerialParam())

saveRDS(bismarkBSseq,file="/home/tander10/sternerlab/brain_methylation/bismarkBSseq_hippocampus.rds")


#2. Script for creating windows and extracting coverage information from our bsseq object

#Modified method for constructing windows
# Load required packages
library(Biostrings)
library(GenomicRanges)

# Assuming you have already loaded the rhesus macaque genome as a DNAStringSet
# Replace 'your_genome' with the actual variable name of your DNAStringSet
genome <- mmul10_fa

# Get the names of all sequences in the genome
seqnames <- names(genome)

# Initialize an empty GRanges object to store the sliding windows
sliding_windows <- GRanges()

# Define the window size
window_size <- 1000

# Iterate through each sequence (chromosome) in the genome
for (seqname in seqnames) {
  # Get the DNAString object for the current sequence
  seq <- genome[[seqname]]

  # Calculate the sliding windows for the current sequence as IRanges
  windows <- IRanges(
    start = seq(1, to = nchar(seq) - window_size + 1, by = window_size),
    end = seq(window_size, to = nchar(seq), by = window_size)
  )

  # Create a GRanges object with the current sequence name
  seq_windows <- GRanges(seqname, ranges = windows)

  # Append the current sequence's sliding windows to the overall sliding_windows object
  sliding_windows <- append(sliding_windows, seq_windows)
}

# Print the first few sliding windows to check
head(sliding_windows)

#Preliminary filtering before extracting coverage information
library(comethyl)
bs <- filterCpGs(bismarkBSseq, cov = 1, perSample = 0.6, file = "/home/tander10/sternerlab/brain_methylation/Filtered_BSseq_hippocampus.rds")

# Calculate the coverage for the windows using getCoverage
coverage <- getCoverage(bs, regions = sliding_windows, type = "Cov", what = "perRegionTotal")

#Get M values for same windows
M_vals <- getCoverage(bs, regions = sliding_windows, type = "M", what = "perRegionTotal")

#Get average methylation estimates for each region
#proportion_meth <- getMeth(bs, regions = sliding_windows, type = "raw", what = "perRegion")

#Generate matrix which contains proportions of methylation at covered sites M/Cov
perc_complete <- M_vals/coverage

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

coverage <- set_granges_as_row_names(coverage, sliding_windows)
M_vals <- set_granges_as_row_names(M_vals, sliding_windows)
#proportion_meth <- set_granges_as_row_names(proportion_meth, sliding_windows)
perc_complete <- set_granges_as_row_names(perc_complete, sliding_windows)

saveRDS(coverage,file="/home/tander10/sternerlab/brain_methylation/RRBScoverage.rds")
saveRDS(M_vals,file="/home/tander10/sternerlab/brain_methylation/RRBSMvals.rds")
#saveRDS(proportion_meth,file="/scratch/tander91/rrbs/RRBSproportion_meth.rds")
saveRDS(perc_complete,file="/home/tander10/sternerlab/brain_methylation/RRBSperc_complete.rds")

#3. Script for filtering out low coverage windows

## Filter for >= 10X median coverage, between 10 and 90% median percent methylation
# Function to filter regions based on coverage and remove NAs
filter_regions <- function(matrix_data, min_coverage) {
  # Calculate the median coverage for each region
  median_coverage <- apply(matrix_data, 1, median, na.rm = TRUE)

  # Filter out regions with less than min_coverage times the median coverage
  filtered_matrix <- matrix_data[median_coverage >= min_coverage, ]

  # Remove NAs from the filtered matrix
  filtered_matrix <- filtered_matrix[complete.cases(filtered_matrix), ]

  return(filtered_matrix)
}

# Filter regions based on coverage
coverage_filtered <- filter_regions(coverage, 10)
meth_filtered <- M_vals[rownames(M_vals) %in% rownames(coverage_filtered), ]
perc_filtered <- perc_complete[rownames(perc_complete) %in% rownames(coverage_filtered), ]
#prop_meth_filtered <- proportion_meth[rownames(proportion_meth) %in% rownames(coverage_filtered), ]

dim(perc_filtered)

filter_methylation <- function(methylation_matrix) {
  # Check if the matrix has row names
  if (is.null(row.names(methylation_matrix))) {
    stop("The input matrix does not have row names.")
  }
  
  # Calculate the mean methylation status for each row
  row_means <- rowMeans(methylation_matrix, na.rm = TRUE)
  
  # Identify rows where the mean methylation status is within the desired range
  keep_rows <- row_means >= 0.1 & row_means <= 0.9
  
  # Filter the matrix based on the identified rows
  filtered_matrix <- methylation_matrix[keep_rows, , drop = FALSE]
  
  # Preserve row names
  row.names(filtered_matrix) <- row.names(methylation_matrix)[keep_rows]
  
  return(filtered_matrix)
}

perc_filtered <- filter_methylation(perc_filtered)
perc_filtered[is.na(perc_filtered)] <- 0
dim(perc_filtered)

saveRDS(coverage_filtered,file="/home/tander10/sternerlab/brain_methylation/RRBScoverage_filtered.rds")
saveRDS(meth_filtered,file="/home/tander10/sternerlab/brain_methylation/RRBSmeth_filtered.rds")
saveRDS(perc_filtered,file="/home/tander10/sternerlab/brain_methylation/RRBSperc_filtered.rds")

#4. Tune alpha: uses the caret package to tune the value of alpha before generating the final clock model.


library(caret)
library(glmnet)

liverMeta <- read.csv(file="/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata.csv")

myControl <- trainControl(method="cv", number = 10, verboseIter = TRUE, search = "grid")
grid <- expand.grid(.alpha=seq(0.1,1,length=5), .lambda=seq(0.001,5,by=0.01))

build_model <- function(x,y){model <- train(x = x, y = y, method = "glmnet", tuneGrid = grid, trControl = myControl)}

model <- build_model(data.matrix(t(perc_filtered)),liverMeta$Age)
plot(model)

unique_rows <- !duplicated(row.names(perc_filtered))
perc_filtered <- perc_filtered[unique_rows, ]

dim(perc_filtered)


# LOOCV: Build Final Model w Leave one out Cross Validation
library(glmnet)

coefs <- data.frame(matrix(nrow = nrow(perc_filtered) + 1, ncol = ncol(perc_filtered)))
row.names(coefs) <- c("intercept", row.names(perc_filtered))

# Create empty lists to store the results
results_predictions <- list()
results_betas <- list()

for (i in 1:ncol(perc_filtered)) {
  train <- perc_filtered[, -i]
  test <- perc_filtered[, i]
  
  trainage <- liverMeta$Age[-i]
  testage <- liverMeta$Age[i]
  testid <- liverMeta$Sample.Name[i]
  trainid <- liverMeta$Sample.Name[-i]
  
  model <- cv.glmnet(data.matrix(t(train)), trainage, nfolds = 10, alpha = 0.325, family = "gaussian")
  
  predicted <- predict(model, newx = t(test), s = model$lambda.min)
  betas <- coef(model, s = model$lambda.min)
  index_beta <- which(!betas[, 1] == 0)
  betacoefs <- betas[index_beta]
  features <- row.names(coefs)[index_beta]
  
  write(rbind(features, betacoefs), ncolumns = 2, file = paste0(testid, "_", i, ".betaCoefficients.coefs"))
  write(c(i, testid, testage, predicted), ncol = 4, file = paste0(testid, "_", i, ".predicted.age"))
  
  # Store the results in the respective lists
  results_predictions[[i]] <- c(i, testid, testage, predicted)
  results_betas[[i]] <- cbind(features, betacoefs)
}

# Write all the predictions to a single file
write.table(do.call(rbind, results_predictions), file = "/home/tander10/sternerlab/brain_methylation/newall_predictions.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write all the beta coefficients to a single file
write.table(do.call(rbind, results_betas), file = "/home/tander10/sternerlab/brain_methylation/all_betaCoefficients.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

predicted_ages <- read.table(file = "/home/tander10/sternerlab/brain_methylation/newall_predictions.txt")

saveRDS(predicted_ages,file="/home/tander10/sternerlab/brain_methylation/hippocampus_predicted_ages.rds")

summary(lm(formula = predicted_ages$V3 ~ predicted_ages$V2))

cor.test(predicted_ages$V3,predicted_ages$V2)

model <- lm(V3 ~ V2, data = predicted_ages)
slope <- coef(model)[2]  # Extract the slope from the model coefficients

# Create the plot
ggplot(predicted_ages, aes(x = V2, y = V3)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "dark green") +
  geom_text(x = max(predicted_ages$V2), y = max(predicted_ages$V3), 
            label = paste("Slope:", round(slope, 2)), hjust = 1, vjust = 1) +
  theme_light()

my_glm_model <- lm(formula = predicted_ages$V3 ~ predicted_ages$V2)

# Step 1: Get the residuals from the GLM model
residuals <- residuals(my_glm_model)

# Step 3: Calculate the median of the absolute residuals (MAD)
mad <- mad(residuals)

# Print the result
print(mad)
