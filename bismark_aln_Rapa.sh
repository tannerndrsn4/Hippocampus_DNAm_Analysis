#!/bin/bash
#SBATCH --partition=memorylong
#SBATCH --account=sternerlab
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --time=10-00:00:00
#SBATCH --array=1-12%12
#SBATCH --job-name=bismark_align_rapamycin
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/Rapamycin/bismark_align_%A_%a.out
#SBATCH --error=/home/tander10/sternerlab/Rapamycin/bismark_align_%A_%a.err
module load miniconda
conda activate bismark

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleList` ## This is how slurm knows to iterate over your sample list; each "array task ID" corresponds to the processing of one file
bismark --genome_folder /home/tander10/sternerlab/Rapamycin/rhesus_genome -1 /home/tander10/sternerlab/Rapamycin/Trim_Galore/${sampleID}*val_1.fq.gz -2 /home/tander10/sternerlab/Rapamycin/Trim_Galore/${sampleID}*val_2.fq.gz -o /home/tander10/sternerlab/Rapamycin/bismark_aligned
conda deactivate
