# AMRIWA
## AMRIWA metagenome analysis

### Analysis using the Snakemake workflow
#### The Snakemake file and the environment `.yml` files are in the `workflow` folder.

#### Create virtual environment
```
conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```
#### Dry-run
```
snakemake --use-conda -np
```
#### Run in batch job
```
snakemake --use-conda -j $SLURM_CPUS_PER_TASK
```

### Analysis run outside the Snakemake workflow are in the `bioinformatics` folder

### The R analysis scripts and related data are in the `RFiles` folder
