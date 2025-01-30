#!/bin/sh
#SBATCH --partition=compute
#SBATCH --job-name=bismark_extract
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/Rapamycin/Bismark_extract_%A_%a.out
#SBATCH --error=/home/tander10/sternerlab/Rapamycin/Bismark_extract_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --account=sternerlab
#SBATCH --array=1-12%12

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleList` ## This is how slurm knows to iterate over your sample list; each "array task ID" corresponds to the processing of one file

module load miniconda
conda activate bismark

bismark_methylation_extractor -p --comprehensive --merge_non_CpG --bedGraph --cytosine_report --gzip --genome_folder /home/tander10/sternerlab/Rapamycin/rhesus_genome -o /home/tander10/sternerlab/Rapamycin/methyl_extract /home/tander10/sternerlab/Rapamycin/bismark_aligned/${sampleID}*.bam

conda deactivate
