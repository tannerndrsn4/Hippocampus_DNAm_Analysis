# Hippocampus_DNAm_Analysis
Differential methylation analysis of RRBS data collected from 96 biobanked rhesus macaque hippocampus samples. 

Data Preprocessing

1. FastQC
2. Trimming w TrimGalore
3. FastQC
4. Genome Index
5. Mapping w Bismark
6. Methylation Extraction

DNA Methylation Analysis

1. Model 1 ~ CpG-based with age as a continuous variable while including sex as a covariate
2. Model 2 ~ DMR-based with age as a continuous variable while including sex as a covariate
3. Model 3 ~ CpG-based ARIMA approach to identify nonlinear patterns in methylation with age while controlling for sex
