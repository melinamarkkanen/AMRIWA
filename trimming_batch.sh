#!/bin/bash
#SBATCH --job-name=Snakemake
#SBATCH --account=Project_2002265  # this is AMRIWA 
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 32
#SBATCH --partition=small

conda activate snakemake
snakemake --use-conda -j $SLURM_CPUS_PER_TASK
