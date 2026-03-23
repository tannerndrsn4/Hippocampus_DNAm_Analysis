# Load required libraries
library(Biostrings)
library(bsseq)
library(ggplot2)
library(parallel)
library(comethyl)
library(limma)

# Load genome reference
#mmul10_fa <- readRDS(file="/projects/sternerlab/shared/Coverage_files/mmul10_fa.rds")
#names(mmul10_fa)[1:21] -> chr
#mmul10_fa <- mmul10_fa[chr]
#loci <- findLoci("CG", subject=mmul10_fa)

# Load Bismark data
#bismarkBSseq <- read.bismark(
#  files = list.files("/home/tander10/sternerlab/brain_methylation/Bismark_output",
#                     full.names = TRUE, pattern = ".cov.gz"),
#  strandCollapse = FALSE,
#  verbose = TRUE,
#  BACKEND = "HDF5Array",
#  loci = loci,
#  rmZeroCov = FALSE,
#  dir = "/projects/sternerlab/shared/Coverage_files/hdf5",
#  replace = TRUE,
#  BPPARAM = BiocParallel::SerialParam()
#)

# Filter CpGs
#bs <- filterCpGs(bismarkBSseq, cov = 5, perSample = 0.6,
#                 file = "/home/tander10/sternerlab/brain_methylation/Filtered_BSseq_hippocampus.rds")
library(HDF5Array)
#bs_file <- "/home/tander10/sternerlab/brain_methylation/results/Filtered_BSseq_Hippocampus_5X.rds"
bs_file <- "/home/tander10/sternerlab/brain_methylation/results/bs_filtered_86samples_hdf5"
#bs <- readRDS(bs_file)
bs <- loadHDF5SummarizedExperiment(bs_file)

# Extract methylation matrix
proportion_meth <- as.matrix(getMeth(bs, type = "raw"))
rownames(proportion_meth) <- paste0(as.character(seqnames(bs)), ":", start(bs))

# Load updated metadata with SVs
# Load updated metadata with SVs
#metadata <- readRDS("/home/tander10/sternerlab/brain_methylation/results/metadata_with_SV1-3.rds")
metadata <- readRDS("/home/tander10/sternerlab/brain_methylation/results/metadata_filtered_86samples_with_SV1-3.rds")
# ✅ Make only 'Age', 'Sex', and 'SubjectID' lowercase
colnames(metadata)[colnames(metadata) == "Age"] <- "age"
colnames(metadata)[colnames(metadata) == "Sex"] <- "sex"
#colnames(metadata)[colnames(metadata) == "SubjectID"] <- "subjectid"

# ✅ Rename 'subjectid' to 'sampid'
#colnames(metadata)[colnames(metadata) == "subjectid"] <- "sampid"
# ✅ Convert sex to a factor
metadata$sex <- as.factor(metadata$sex)
colnames(proportion_meth) <- trimws(as.character(metadata$sampid))

# Remove rows with NA values
proportion_meth_clean <- na.omit(proportion_meth)

# Match samples
matching_samples <- intersect(colnames(proportion_meth_clean), metadata$sampid)
proportion_meth_clean <- proportion_meth_clean[, matching_samples, drop = FALSE]
metadata <- metadata[metadata$sampid %in% matching_samples, , drop = FALSE]

# Adjust for Sex + SV1–3
design <- model.matrix(~ sex + SV1 + SV2 + SV3, data = metadata)
residual_meth <- removeBatchEffect(proportion_meth_clean, design = design)
dim(residual_meth)

# Save residuals
saveRDS(residual_meth, file = "Percent_Methylated_Residuals_SexSVAdj_Limma.rds")
saveRDS(proportion_meth_clean, file = "/home/tander10/sternerlab/brain_methylation/Percent_Methylated_hippocampus.rds")
write.csv(metadata, file = "/home/tander10/sternerlab/brain_methylation/Brain_RRBS_metadata_modified.csv", row.names = FALSE)

# Define ARIMA function
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
# (Assume arima.selection is defined above in full)

# Run ARIMA on CpG-level
#aging_arima_cpg <- arima.selection(adj.data = residual_meth, metadata = metadata)
#save(aging_arima_cpg, file = "/home/tander10/sternerlab/brain_methylation/Hippocampus_arima_CpG_level_SVAdj.RData")

### ---------------- REGION-LEVEL ARIMA ----------------

# Region identification
meth_cov <- getCoverage(bs, type = "Cov")
CG_positions <- start(rowRanges(bs))
CG_chr <- as.character(seqnames(rowRanges(bs)))
dummy_stat <- rep(1, length(CG_positions))
region_data <- bsseq:::regionFinder3(
  x = dummy_stat,
  chr = CG_chr,
  positions = CG_positions,
  maxGap = 1000,
  verbose = TRUE
)[["up"]]

region_data$RegionID <- paste0(region_data$chr, "_", region_data$start, "_", region_data$end)

# --- NEW: add length and CpG counts ---
library(GenomicRanges)
region_data$length <- region_data$end - region_data$start + 1

cpg_gr <- GRanges(seqnames = CG_chr, ranges = IRanges(start = CG_positions, width = 1))
region_gr <- GRanges(seqnames = region_data$chr,
                     ranges = IRanges(start = region_data$start, end = region_data$end))
region_data$num_CpGs <- countOverlaps(region_gr, cpg_gr)

# --- NEW: filter out regions with <2 CpGs ---
region_data <- region_data[region_data$num_CpGs >= 3, ]

# --- NEW: summary statistics ---
cat("CpGs per region:\n")
print(summary(region_data$num_CpGs))
cat("Min CpGs:", min(region_data$num_CpGs), "\n")
cat("Max CpGs:", max(region_data$num_CpGs), "\n\n")

cat("Region lengths (bp):\n")
print(summary(region_data$length))
cat("Min length:", min(region_data$length), "\n")
cat("Max length:", max(region_data$length), "\n")

# Extract region-level methylation
coverage_region <- getCoverage(bs, type = "Cov", regions = region_data, what = "perRegionTotal")
methylation_region <- getCoverage(bs, type = "M", regions = region_data, what = "perRegionTotal")
perc_meth_region <- methylation_region / coverage_region

# Match column names to sample IDs from metadata
colnames(perc_meth_region) <- trimws(as.character(metadata$sampid))
rownames(perc_meth_region) <- region_data$RegionID

# Filter regions by coverage and variability
filter_regions <- function(matrix_data, min_coverage) {
  med_cov <- apply(matrix_data, 1, median, na.rm = TRUE)
  matrix_data[med_cov >= min_coverage, ]
}
coverage_filtered <- filter_regions(coverage_region, 0)

perc_filtered <- perc_meth_region
perc_filtered <- perc_filtered[
  rowMeans(perc_filtered, na.rm = TRUE) >= 0.1 &
  rowMeans(perc_filtered, na.rm = TRUE) <= 0.9, , drop = FALSE]
perc_filtered <- perc_filtered[complete.cases(perc_filtered), ]
dim(perc_filtered)

# Match metadata to region-level matrix
matching_samples_region <- intersect(colnames(perc_filtered), metadata$sampid)
perc_filtered <- perc_filtered[, matching_samples_region, drop = FALSE]
metadata_region <- metadata[metadata$sampid %in% matching_samples_region, , drop = FALSE]
dim(metadata_region)

# Adjust region-level methylation for Sex + SV1–3
residual_region <- removeBatchEffect(
  perc_filtered,
  design = model.matrix(~ sex + SV1 + SV2 + SV3, data = metadata_region)
)

# Save region-level residuals
saveRDS(residual_region,
        "/home/tander10/sternerlab/brain_methylation/Percent_Methylated_Regions_Residuals_SexSVAdj_Limma.rds")

# Run ARIMA on regions
aging_arima_region <- arima.selection(adj.data = residual_region, metadata = metadata_region)
save(aging_arima_region,
     file = "/home/tander10/sternerlab/brain_methylation/Hippocampus_arima_Region_level_SVAdj.RData")

# ============================================================
# SEX-SPECIFIC REGION-LEVEL ARIMA
# ============================================================

output_dir_arima_sex <- "/home/tander10/sternerlab/brain_methylation/results_arima_region_sex"
dir.create(output_dir_arima_sex, showWarnings = FALSE, recursive = TRUE)

# Ensure expected column names exist (your script already makes sex a factor)
stopifnot(all(c("sampid", "age", "sex") %in% colnames(metadata_region)))

# Optional: ensure consistent sample ordering between matrix and metadata
metadata_region <- metadata_region[match(colnames(perc_filtered), metadata_region$sampid), ]
stopifnot(all(colnames(perc_filtered) == metadata_region$sampid))

# Split sample IDs by sex
metadata_region$sampid <- as.character(metadata_region$sampid)
colnames(perc_filtered) <- as.character(colnames(perc_filtered))

sex_groups <- split(metadata_region$sampid, metadata_region$sex)
head(sex_groups)
print(lapply(sex_groups, length))

head(perc_filtered)

# Loop over sexes
for (sx in names(sex_groups)) {
  cat("\n==============================\n")
  cat(" Running REGION-level ARIMA for sex:", sx, "\n")
  cat("==============================\n")

  sx_dir <- file.path(output_dir_arima_sex, paste0("ARIMA_Region_", sx))
  dir.create(sx_dir, showWarnings = FALSE, recursive = TRUE)

  sampids_sx <- sex_groups[[sx]]

  # Subset region methylation matrix + metadata
  perc_sx <- perc_filtered[, sampids_sx, drop = FALSE]
  meta_sx <- metadata_region[metadata_region$sampid %in% sampids_sx, , drop = FALSE]

  # Re-order metadata to match matrix columns exactly
  meta_sx <- meta_sx[match(colnames(perc_sx), meta_sx$sampid), , drop = FALSE]
  stopifnot(all(colnames(perc_sx) == meta_sx$sampid))

  # Adjust within sex: SV1–SV3 only (NO sex term)
  design_sx <- model.matrix(~ SV1 + SV2 + SV3, data = meta_sx)

  residual_region_sx <- removeBatchEffect(perc_sx, design = design_sx)

  # Save residuals
  saveRDS(
    residual_region_sx,
    file = file.path(sx_dir, paste0("Percent_Methylated_Regions_Residuals_SVAdj_", sx, ".rds"))
  )

  # Run ARIMA (same function as before)
  aging_arima_region_sx <- arima.selection(
    adj.data = residual_region_sx,
    metadata = meta_sx,
    dataType = "rna",   # ok to leave; only affects printed labels in your function
    n.cores = 20
  )

  # Save ARIMA results
  save(
    aging_arima_region_sx,
    file = file.path(sx_dir, paste0("Hippocampus_arima_Region_level_", sx, "_SVAdj.RData"))
  )

  cat("✅ Completed REGION-level ARIMA for", sx, "\n")
}

cat("\n🎉 Sex-specific REGION-level ARIMA complete!\n")
