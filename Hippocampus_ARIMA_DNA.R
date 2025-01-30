# 1. Script for loading in your data and creating a bsseq object

library(Biostrings)
library(bsseq)
library(ggplot2)
library(parallel)
library(comethyl)

# Load genome reference
mmul10_fa <- readRDS(file="/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")

# Filter for the first 20 chromosomes (chr1-chrX)
names(mmul10_fa)[1:21] -> chr
mmul10_fa <- mmul10_fa[chr]

# Find CG sites in the genome
loci <- findLoci("CG", subject=mmul10_fa)

# Read Bismark data
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

# Preliminary filtering before extracting coverage information
bs <- filterCpGs(bismarkBSseq, cov = 10, perSample = 0.6, 
                 file = "/home/tander10/sternerlab/brain_methylation/Filtered_BSseq_hippocampus.rds")

# Get average methylation estimates for each region
proportion_meth <- as.matrix(getMeth(bs, type = "raw"))

# Debugging: Print initial dimensions
cat("Initial dimensions of proportion_meth:", dim(proportion_meth), "\n")

# **Assign row names based on CpG coordinates**
cat("Assigning row names to proportion_meth based on genomic coordinates...\n")

# Extract chromosome and position from the BSseq object
chr <- as.character(seqnames(bs))   # Extract chromosome names
pos <- start(bs)                    # Extract CpG start positions

# Construct row names in the format "chr:pos" (e.g., "chr1:12345")
rownames(proportion_meth) <- paste0(chr, ":", pos)

# Debugging: Check first few row names
cat("First few row names in proportion_meth:\n")
print(head(rownames(proportion_meth)))

# Read metadata file
metadata <- read.csv(file='/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata.csv')

# Make all column names lowercase
colnames(metadata) <- tolower(colnames(metadata))

# Rename "subjectid" to "sampid"
colnames(metadata)[colnames(metadata) == "subjectid"] <- "sampid"

# Ensure sampid has no leading/trailing spaces
metadata$sampid <- trimws(as.character(metadata$sampid))

# Assign sample IDs to proportion_meth columns
colnames(proportion_meth) <- trimws(as.character(metadata$sampid))

# Remove rows with NA values in the matrix
cat("Checking for NA values in proportion_meth...\n")
print(sum(is.na(proportion_meth)))  # Total NA count
if (anyNA(proportion_meth)) {
  cat("⚠️ Warning: NA values detected! Removing rows with NAs...\n")
  proportion_meth_clean <- na.omit(proportion_meth)  # Drop rows with any NA
} else {
  proportion_meth_clean <- proportion_meth
}

# Keep only matching samples in matrix and metadata
matching_samples <- intersect(colnames(proportion_meth_clean), metadata$sampid)
proportion_meth_clean <- proportion_meth_clean[, matching_samples, drop = FALSE]
metadata <- metadata[metadata$sampid %in% matching_samples, , drop = FALSE]

# Debugging: Final dimensions check
cat("Final dimensions of proportion_meth_clean: ", dim(proportion_meth_clean), "\n")
cat("Final metadata rows: ", nrow(metadata), "\n")

# Save cleaned data
saveRDS(proportion_meth_clean, file = "/home/tander10/sternerlab/brain_methylation/Percent_Methylated_hippocampus.rds")
write.csv(metadata, file="/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata_modified.csv", row.names = FALSE)

# ARIMA fxn - from Marquez 2020
arima.selection <- function(adj.data=NULL,
                            metadata=NULL,
                            evalage=NULL,
                            use.boxcox=FALSE, # include selected.sex=NULL if you want to run ARIMA by sexes individually 
                            use.portmanteau=TRUE,
                            loess.span=0.75,
                            dataType=c("rna"),
                            n.cores=20) {
  require(zoo)
  require(forecast)
  require(reshape2)
  require(parallel)
  require(dplyr)
  
  age=with(metadata,setNames(age,sampid)) # pairs LID with age 
  sex=with(metadata,setNames(sex,sampid)) # pairs LID with sex (but this code does not separate GE by sex, original Marquez code does)
  
  cat("Building time series object\n")
  age.sexsel <- age[names(sex)]
  age.sexsel.sorted <- age.sexsel[order(age.sexsel)] # orders by age
  adj.data.sexsel <- adj.data[,names(sex)]
  
  selpeaks.sexsel <- t(adj.data.sexsel)
  age.counts <- as.data.frame(table(age.sexsel),stringsAsFactors = F)
  selpeaks.sexsel.unique <- do.call(rbind,lapply(setNames(as.list(seq_len(nrow(age.counts))),age.counts$age.sexsel),function(n) {
    if (age.counts[n,"Freq"]==1) {
      estage <- selpeaks.sexsel[age.sexsel==age.counts[n,"age.sexsel"],]
    } else {
      estage <- colMeans(selpeaks.sexsel[age.sexsel==age.counts[n,"age.sexsel"],])
    }
    return(estage)}))
  age.sexsel.unique <- as.numeric(age.counts$age.sexsel)
  
  # Create zoo/ts objects, average data for same ages, interpolate unsampled ages
  zoo.sexsel <- zoo(selpeaks.sexsel.unique,order.by=age.sexsel.unique,frequency = 1)
  ts.sexsel <- mclapply(zoo.sexsel,ts,mc.cores = n.cores)
  # plot.ts(ts.sexsel[[1]],type="b") # Visualization of a single TS
  if (use.boxcox) {
    cat("Calculating Box-Cox Lambda using the log-likelihood method\n")
    bclambda.sexsel <- sapply(ts.sexsel,BoxCox.lambda,method="loglik",upper=10)
  } else {
    bclambda.sexsel = NULL
  }
  gc()
  
  cat("Fitting auto-arima models on",n.cores,"CPU cores\n")
  
  # Builds an ARIMA model to test for trends in each time series, then selects those for which the trend is non-negligible (stationarity cannot be rejected AND no AR or MA terms are chosen)
  arima.sexsel <- mclapply(setNames(seq_len(length(ts.sexsel)),names(ts.sexsel)),function(n) auto.arima(ts.sexsel[[n]],seasonal = F,allowmean = T,allowdrift = T,lambda = bclambda.sexsel), mc.cores = n.cores)
  gc()
  if (use.portmanteau) {
    # portmanteau tests on RESIDUALS (Ljung-Box) - models with autocorrelated residuals are removed as well
    cat("Computing portmanteau tests on ARIMA model residuals using the Ljung-Box algorithm, on",n.cores,"CPU cores\n")
    portmanteau.sexsel <- unlist(mclapply(arima.sexsel,function(x) {
      df=sum(x$arma[1:2])
      pmin <- Box.test(x$residuals,lag=20,type="Ljung-Box",fitdf = df)$p.value
      return(pmin)
    },mc.cores = n.cores))
    cat("Proceeding to filter out",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"with stationary models and autocorrelated residuals\n")
    arima.sexsel.select <- unlist(lapply(arima.sexsel,function(x) sum(x$arma[1:2])>0 & x$arma[6]>0)) & portmanteau.sexsel>0.05
  } else {
    cat("Proceeding to filter out",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"with stationary models and autocorrelated residuals\n")
    arima.sexsel.select <- unlist(lapply(arima.sexsel,function(x) sum(x$arma[1:2])>0 & x$arma[6]>0))
    portmanteau.sexsel = NULL
  }
  gc()
  
  arima.sexsel.nonzero <- arima.sexsel[arima.sexsel.select]
  ts.sexsel.nonzero <- ts.sexsel[arima.sexsel.select]
  fitted.sexsel.nonzero <- t(as.data.frame(sapply(arima.sexsel.nonzero,`[[`,"fitted"),row.names = as.integer(attr(arima.sexsel.nonzero[[1]]$x,"index"))))
  cat("Out of",length(arima.sexsel),ifelse(dataType[[1]]=="atac","peaks,","transcripts,"),length(ts.sexsel.nonzero),"trendy and well-fitted",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"were selected for further analysis\n")
  gc()
  
  ## Uses local fit (LOESS) on sex-union nonzero trends to predict chromatin values at each age in the span, for sex comparisons
  cat("Computing predicting",ifelse(dataType[[1]]=="atac","peaks","transcripts"),"using LOESS interpolation over ARIMA-fitted data on",n.cores,"cores\n")
  loess.sexsel.nonzero <- mclapply(arima.sexsel.nonzero,function(A) loess(as.vector(A$fitted)~attr(A$x,"index"),span = loess.span,family = "symmetric"),mc.cores = n.cores)
  gc()
  predicted.sexsel.nonzero <- t(as.data.frame(mclapply(loess.sexsel.nonzero,predict,evalage,mc.cores = n.cores),row.names = evalage))
  gc()
  
  output <- list(ts.complete=ts.sexsel,
                 arima.complete=arima.sexsel,
                 arima.select=arima.sexsel.select,
                 loess.nonzero=loess.sexsel.nonzero,
                 fitted.nonzero=fitted.sexsel.nonzero,
                 predicted.nonzero=predicted.sexsel.nonzero,
                 BoxCox_lambda=bclambda.sexsel,
                 portmanteau.pvalues=portmanteau.sexsel)
  
  cat("Done.\n")
  return(output)
}

## Run with: 
aging_arima <- arima.selection(adj.data = resid_counts, metadata = covariates)
save(aging_arima, file = "aging_miRNA_arima.RData")

# Run ARIMA
aging_arima <- arima.selection(adj.data = proportion_meth_clean, metadata = metadata)
save(aging_arima, file = "/home/tander10/sternerlab/brain_methylation/Hippocampus_arima_DNAmethylation.RData")