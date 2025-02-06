#!/bin/sh
#SBATCH --licenses=common
#SBATCH --partition=samodha
#SBATCH --ntasks-per-node=4
#SBATCH --mem=90GB
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --job-name=16S_analysis
#SBATCH --error=dada2_16S_analysis.%J.stdout
#SBATCH --output=dada2_16S_analysis.%J.stderr
#SBATCH --mail-user=adamyazori@gmail.com
#SBATCH --mail-type=ALL


#module load R/4.1
#module load mothur
conda activate R

#analysis
Rscript 16S.R

#chmod +x mothur.sh

#bash mothur.sh
