# AMRIWA
AMRIWA metagenome analysis

# Using the Snakemake workflow

Make a snakemake virtual environment using Mamba.  
First install Mamba and then Snakemake.
```
conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```

The Snakemake file and the virtual environment `.yml` files are in the `workflow` folder.  
The virtual environments for QC and cutadapt can be created first with:
```
snakemake --use-conda --conda-create-envs-only --cores 1
```

Then check with dry-run that everything works.
```
snakemake --use-conda -np
```

If everything is ok, proceed to trimming the reads. All data should be in the `data` folder (which is not included in git).  
And after trimming trimmed data will be in `trimmed_data` folder (also not included in git).   
FastQC and MultiQC files for raw and trimmed data can be found from corresponding folders as well.
```
snakemake --use-conda -j 32
```

Probably easiest to run as a batch job in Puhti.
Example bacth file:
```
#!/bin/bash
#SBATCH --job-name=Snakemake
#SBATCH --account=Project_2002265
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 32
#SBATCH --partition=small

conda activate snakemake
snakemake --use-conda -j $SLURM_CPUS_PER_TASK
```

