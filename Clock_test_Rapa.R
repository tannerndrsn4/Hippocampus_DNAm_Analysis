library(glmnet)
library(bsseq)
library(GenomicRanges)


# # 1. Load in ONPRC clock data 
# 
perc_filtered <- readRDS(file="/projects/sternerlab/shared/Coverage_files/ONPRC_perc_filtered.rds")
perc_filtered[is.na(perc_filtered)] <- 0
liverMeta <- readRDS(file="/projects/sternerlab/shared/Coverage_files/ONPRCmetadata.rds")

# # 2. Construct your bsseq object with test data
# 
TestMeta <- readRDS(file="/projects/sternerlab/shared/Coverage_files/CRMetadata_20211228.rds")
#Fasta of the rhesus macaque reference genome
mmul10_fa <- readRDS(file="/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")

## filter for the first 21 chromosomes (chr1-chrX)
names(mmul10_fa)[1:21]->chr
mmul10_fa=mmul10_fa[chr]

## find CG sites in the genome
loci = findLoci("CG", subject=mmul10_fa)

bismarkBSseq_test <- read.bismark(files = list.files("/projects/sternerlab/shared/Coverage_files/Wisc_extract_cov",
                                                full.names = T, pattern=".cov.gz"),
                             strandCollapse = FALSE,
                             verbose = TRUE,
                             BACKEND = "HDF5Array",
                             loci = loci,
                             rmZeroCov = FALSE,
                             dir="/projects/sternerlab/shared/Coverage_files/hdf5",
                             replace = TRUE,
                             BPPARAM = BiocParallel::SerialParam())

saveRDS(bismarkBSseq_test,file="/projects/sternerlab/shared/Coverage_files/bismarkBSseq_Wisconsin.rds")

# 3. Extract coverage information for the windows of interest in our test data

# Extract chromosome, start, and end positions from row names
coordinates <- strsplit(rownames(perc_filtered), "_")
chromosomes <- sapply(coordinates, function(coord) coord[1])  # Keep chromosome names as characters
starts <- sapply(coordinates, function(coord) as.numeric(coord[2]))
ends <- sapply(coordinates, function(coord) as.numeric(coord[3]))

# Create the GRanges object
reg <- GRanges(seqnames = chromosomes, ranges = IRanges(start = starts, end = ends))

coverage_test <- getCoverage(bismarkBSseq_test, regions = reg, type = "Cov", what = "perRegionTotal")
M_vals_test <- getCoverage(bismarkBSseq_test, regions = reg, type = "M", what = "perRegionTotal")
#Generate matrix which contains proportions of methylation at covered sites M/Cov
perc_complete_test <- M_vals_test/coverage_test

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

# 4. Use filter function to remove any regions for which there are missing data

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
coverage_filtered_test <- filter_regions(coverage_test, 0)
meth_filtered_test <- M_vals_test[rownames(M_vals_test) %in% rownames(coverage_filtered_test), ]
perc_filtered_test <- perc_complete_test[rownames(perc_complete_test) %in% rownames(coverage_filtered_test), ]

# saveRDS(perc_filtered_test,file="/projects/sternerlab/shared/Coverage_files/perc_filtered_Wisconsin.rds")

#perc_filtered_test <- readRDS(file="/projects/sternerlab/shared/Coverage_files/perc_filtered_Wisconsin.rds")
perc_filtered_test[is.na(perc_filtered_test)] <- 0

# 5. Now to run our test data through the clock -- there are two methods for this, the first being the cor.test and the second using the predict function 

norm_train <- apply(perc_filtered, 1, function(x) {return(qqnorm(x, plot = FALSE)$x)})
norm_train <- t(norm_train)
norm_train <- as.matrix(norm_train)

norm_train2 <- apply(norm_train, 2, function(x) {return(qqnorm(x, plot = FALSE)$x)})
norm_train2 <- t(norm_train2)
colnames(norm_train2) <- colnames(norm_train)


norm_test<-apply(perc_filtered_test,1,function(x){return(qqnorm(x,plot=F)$x)})
norm_test <- t(norm_test)
norm_test <- as.matrix(norm_test)

norm_test2<-as.matrix(apply(norm_test,2,function(x){return(qqnorm(x,plot=F)$x)}))
norm_test2<-t(norm_test2)
colnames(norm_test2) <- colnames(norm_test)


# #Build final model: extracts the final window position identifiers and beta coefficients for features in the clock and saves them as a txt file.
coefs <- data.frame(matrix(nrow=nrow(perc_filtered)+1,ncol=ncol(perc_filtered)))
rownames(coefs) <- c("intercept",rownames(perc_filtered))
set.seed(9)
finalModel <- cv.glmnet(data.matrix(t(perc_filtered)),y=liverMeta$Age,family = "gaussian",alpha = 0.1)


## get non-zero coefficients from model
betas <- coef(finalModel)
betasWithFeatures <- data.frame(
  features = betas@Dimnames[[1]][ which(betas != 0 ) ], #intercept included
  coefs    = betas              [ which(betas != 0 ) ]  #intercept included
)


## save coefficients and features to be used for the clock
write.table(betasWithFeatures,file="/projects/sternerlab/shared/Coverage_files/ClockBetas_Wisconsin.txt",quote = F,sep="\t",row.names = F,col.names = F)

#library(tidyr)
#betasWithFeatures <- read.delim("/projects/sternerlab/shared/Coverage_files/10003_13.betaCoefficients.coefs", sep="\t", header=FALSE)
#betasWithFeatures <- separate(betasWithFeatures, V1, into = c("col1", "col2"), sep = " ")

#colnames(betasWithFeatures) <- c("features", "coefs")
#betasWithFeatures$coefs <- as.numeric(betasWithFeatures$coefs)
#head(betasWithFeatures)

# #6. Predict R:
# ## save intercept value (optional; can just reference original text file)
intercept <- betasWithFeatures$coefs[1]
## Remove intercept from feature IDs (chrom_pos)
clock_features <- betasWithFeatures$features[-1]
## remove intercept from beta coefficients
clock_coefs <- betasWithFeatures$coefs[-1]
## subset percent methylation matrix (perc_filtered) to retain only windows in the clock in the perc meth matrix; format of rownames must match format of features ("chr_pos")
reduced_clock <- perc_filtered_test[rownames(perc_filtered_test) %in% clock_features,]

## calculate sum of weighted regression coefficients
combine <- apply(reduced_clock, 2, function(x){x*clock_coefs})
colsums_combine <- as.data.frame(colSums(combine))+intercept ## add intercept to totals
colsums_combine
saveRDS(colsums_combine, file = "/projects/sternerlab/shared/Coverage_files/LOOCV_predictedages.rds")

dim(perc_filtered)
dim(perc_filtered_test)

# #Making the plot
# library(readr)

# predicted_ages <- data.frame(
#   Number = length(TestMeta$Age_at_Death),
#   Sample_ID = TestMeta$Animal_ID,
#   Chronological_Age = TestMeta$Age_at_Death,
#   Predicted_Age = colsums_combine$`colSums(combine)`
#   )
# 
# saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/predicted_ages_Wisconsin.rds")

## examine strength of correlation (should be identical or very similar to glmnet's automatic results)
# cor.test(colsums_combine$`colSums(combine)`,TestMeta$Age_at_Death)

# 6. Other approach using the predict function

#finalModel <- cv.glmnet(data.matrix(t(perc_filtered)),y=liverMeta$Age,family = "gaussian",alpha = 0.55)
finalModel <- cv.glmnet(data.matrix(t(perc_filtered)),y=liverMeta$Age,family = "gaussian",alpha = 0.325)

#If you are missing coverage information for any respective windows in the test dataset, you may need to filter the ONPRC dataset to include only windows that match the test dataset

predicted_ages <- predict(finalModel, t(perc_filtered_test), s = "lambda.min")
#predicted_ages <- predict(finalModel, t(perc_filtered_test))
predicted_ages
saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/predicted_ages_Wisconsin_2.rds")

finalModel <- cv.glmnet(data.matrix(t(norm_train2)),y=liverMeta$Age,family = "gaussian",alpha = 0.325)
predicted_ages <- predict(finalModel, t(norm_test2))
predicted_ages
saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/newpredicted_ages_Wisconsin_norm.325.rds")

# predicted_ages <- as.data.frame(predicted_ages)
#predicted_ages$chrono <- TestMeta$Age.at.death..years.

#saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/predicted_ages_Wisconsin_2.rds")

#lm_model <- lm(lambda.min ~ chrono, data = predicted_ages)
#summary(lm_model)

# Extract the slope from the linear regression model
#slope <- coef(lm_model)[2]  # Extract the coefficient for Chronological_Age

# Create the scatter plot with the linear regression line and slope annotation
#library(ggplot2)

# p <- ggplot(predicted_ages, aes(x = chrono, y = lambda.min)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "dark green", se = FALSE) +  # Add linear regression line
#   geom_text(x = max(predicted_ages$chrono), y = max(predicted_ages$lambda.min),
#             label = paste("Slope = ", round(slope, 2)), hjust = 1, vjust = 1) +  # Add slope annotation
#   theme_light()

#ggsave(filename = "Wisconsin_agepredictions.pdf", p)







# finalModel <- cv.glmnet(data.matrix(t(norm_train2)),y=liverMeta$Age,family = "gaussian",alpha = 0.1)
# predicted_ages <- predict(finalModel, t(norm_test2))
# predicted_ages
# saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/newpredicted_ages_Wisconsin_norm.1.rds")
# 
# finalModel <- cv.glmnet(data.matrix(t(norm_train2)),y=liverMeta$Age,family = "gaussian",alpha = 0.9)
# predicted_ages <- predict(finalModel, t(norm_test2))
# predicted_ages
# saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/newpredicted_ages_Wisconsin_norm.9.rds")
# 
# finalModel <- cv.glmnet(data.matrix(t(norm_train2)),y=liverMeta$Age,family = "gaussian",alpha = 1.0)
# predicted_ages <- predict(finalModel, t(norm_test2))
# predicted_ages
# saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/newpredicted_ages_Wisconsin_norm1.0.rds")
# 
# finalModel <- cv.glmnet(data.matrix(t(norm_train2)),y=liverMeta$Age,family = "gaussian",alpha = 0.5)
# predicted_ages <- predict(finalModel, t(norm_test2))
# predicted_ages
# saveRDS(predicted_ages,file="/projects/sternerlab/shared/Coverage_files/newpredicted_ages_Wisconsin_norm.5.rds")
# 
# #Code for making plots 
# 
# predicted_ages <- as.data.frame(predicted_ages)
# 
# predicted_ages$group <- Wisc_metadat$Group
# predicted_ages$chrono <- Wisc_metadat$Age_at_Death
# 
# # Run linear models within each group
# lm_model_C <- lm(lambda.min ~ chrono, data = subset(predicted_ages, group == "C"))
# lm_model_R <- lm(lambda.min ~ chrono, data = subset(predicted_ages, group == "R"))
# 
# # Extract slopes from the linear regression models
# slope_C <- coef(lm_model_C)[2]  # Extract the coefficient for Chronological_Age in group C
# slope_R <- coef(lm_model_R)[2]  # Extract the coefficient for Chronological_Age in group R
# 
# # Create the scatter plot with separate linear regression lines and slope annotations for each group
# library(ggplot2)
# 
# p <- ggplot(predicted_ages, aes(x = chrono, y = lambda.min, color = group)) +
#   geom_point() +
#   geom_smooth(aes(group = group), method = "lm", se = FALSE) +  # Add separate linear regression lines for each group
#   geom_text(data = data.frame(group = c("C", "R"), slope = c(slope_C, slope_R)),
#             aes(x = Inf, y = Inf, label = paste("Slope = ", round(slope, 2))),
#             hjust = 1, vjust = 1, color = c("dark green", "blue")) +  # Add slope annotations for each group with corresponding colors
#   theme_light()
# 
# # Display the plot
# print(p)
# 
# # Step 1: Get the residuals from the GLM model
# residuals <- residuals(lm_model_C)
# 
# # Step 2: Take the absolute value of the residuals
# abs_residuals <- abs(residuals)
# 
# # Step 3: Calculate the median of the absolute residuals (MAD)
# mad <- median(abs_residuals)
# 
# # Print the result
# print(mad)
