# Hippocampus_DNAm_Analysis
Differential methylation analysis of RRBS data collected from 86 biobanked rhesus macaque hippocampus samples. Includes pipelines for DNAm analyses and code used for figure generation.  

## Data Preprocessing

1. [FastQC](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/FastQC_rapa.sh)
2. [Trimming](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/TrimGalore_rapa.sh) w TrimGalore
3. [FastQC](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/FastQC_trimmed_rapa.sh) on trimmed reads
4. [Genome Index](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/bismark_genome_preparationRapa.sh)
5. [Mapping](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/bismark_aln_Rapa.sh) w Bismark
6. [Methylation Extraction](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/Methylation_extract_Rapa.sh)

## DNA Methylation Analysis

1. CpG-based [differential methylation](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/Hippocampus_DMR_Analyses.R) and [variance](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/Hippocampus_Variance_Analyses.r) analysis
2. Region-based [identification of age-associated regions](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/Hippocampus_ARIMA_DNA.R) using ARIMA modeling and [clustering](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/Hippocampus_ARIMA_Visualization_REGION.r)
3. Genic [integration with matched RNAseq](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/IntersectionwRNAseq_Approach1.R)
4. Intergenic [intersection analyses](https://github.com/tannerndrsn4/Hippocampus_DNAm_Analysis/blob/main/Enhancer_TE_analysis.R) to identify overlap with transposable elements and human enhancers 
