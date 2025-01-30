#!/bin/sh
#SBATCH --partition=compute
#SBATCH --job-name=FastQC_trimmed_rapa
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/Rapamycin/FastQC_%A_%a.out ## %A_%a means add "slurmJobID_ArrayID" to the end of the file name
#SBATCH --error=/home/tander10/sternerlab/Rapamycin/FastQC_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --account=sternerlab
#SBATCH --array=1-12%12

## Above: For the array, start with just --array=1 to make sure it works; it shouldn't take more than an hour at absolute most; then you can run: --array=2-192%100 ## This means run all samples in groups of 100 (so roughly two batches); you can modify as you wish (e.g., if you have 50 files you want to run the array job on, then you would change it to --array=1-50%50)


## Below, 'sampleList' is a separate text file you need to make beforehand, containing the unique names of your files saved in the same directory. The names in 'SampleList' *should not* include the .fastq suffix (otherwise Slurm will think your file is named 'myfile.fastq.fastq' because of the way the place holder 'sampleID' is formatted; if you send me your file names, I can show you how to automatically generate this list using the `cut` command if you're not already familiar; don't do it by hand!


sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleList` ## This is how slurm knows to iterate over your sample list; each "array task ID" corresponds to the processing of one file

module load fastqc/0.11.5

fastqc -o /home/tander10/sternerlab/Rapamycin/Fastqc_trimmed /home/tander10/sternerlab/Rapamycin/Trim_Galore/${sampleID}*R1*fq.gz /home/tander10/sternerlab/Rapamycin/Trim_Galore/${sampleID}*R2*fq.gz
