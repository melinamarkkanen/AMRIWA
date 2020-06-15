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
The virtual environments for QC and cutadaopt can be created first with:
```
snakemake --use-conda --conda--create-envs-only
