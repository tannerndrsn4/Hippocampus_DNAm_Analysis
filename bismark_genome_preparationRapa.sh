#!/bin/sh
#SBATCH --partition=compute
#SBATCH --job-name=bismark_genome_preparation
#SBATCH --mail-user=tander10@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/tander10/sternerlab/Rapamycin/bismark_genome_preparation.output
#SBATCH --error=/home/tander10/sternerlab/Rapamycin/bismark_genome_preparation.err
#SBATCH --time=24:00:00
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --account=sternerlab

module load miniconda
conda activate bismark

bismark_genome_preparation --verbose /home/tander10/sternerlab/Rapamycin/rhesus_genome

conda deactivate
