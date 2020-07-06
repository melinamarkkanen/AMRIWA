# AMRIWA
AMRIWA metagenome analysis

# Using the Snakemake workflow

Make a snakemake virtual environment in Puhti using Mamba.
```
conda create -n mamba
source activate mamba
conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```

The Snakemake file and the virtual environment `.yml` files are in the `workflow` folder.  
The virtual environments can be created first with:
```
snakemake --use-conda --conda-create-envs-only --cores 1
```

Then check with dry-run that everything works.
```
snakemake --use-conda -np
```
Batch job for dry run in Puhti
```
#!/bin/bash
#SBATCH --job-name=Snakemake
#SBATCH --account=Project_2002265
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --partition=test
#SBATCH --output=snakemake_out_%j.txt

module load bioconda/3
source activate mamba
source activate snakemake

cd /scratch/project_2002265/markkan5/AMRIWA/workflow

snakemake --use-conda -j $SLURM_CPUS_PER_TASK -np

conda deactivate
conda deactivate
```

Batch job in Puhti. Add "--latency-wait" if needed.

```
#!/bin/bash
#SBATCH --job-name=Snakemake
#SBATCH --account=Project_2002265
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 32
#SBATCH --partition=small
#SBATCH --output=snakemake_out_%j.txt
#SBATCH --error=snakemake_err_%j.txt

module load bioconda/3
source activate mamba
source activate snakemake

cd /scratch/project_2002265/markkan5/AMRIWA

snakemake --use-conda -j $SLURM_CPUS_PER_TASK

conda deactivate
conda deactivate
```
Download and uncompress data:

First set password for Allas for project. Must be run every time prior to actual batch job.
```
module load allas
allas-conf -k project_2002265
```

Then download and uncompress data with batch job:
```
#!/bin/bash -l
#SBATCH -J download_data
#SBATCH --account=project_2002265
#SBATCH -o download_data_out_%A_%a.txt
#SBATCH -e download_data_err_%A_%a.txt
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --array=1-2
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -p small

# go to target dir
cd /scratch/project_2002265/markkan5/AMRIWA/workflow/data

# set environmental variable for sample names
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2002265/markkan5/AMRIWA/sample_names.txt)

# make sure the connection to Allas is open
source /appl/opt/allas-cli-utils/allas_conf -f -k $OS_PROJECT_NAME

# download the whole bucket
swift download 2002265_Melina_Markkanen_AMRIWA_metagenomes

# remove the files that are not used here
rm *_ameta
rm FASTQC.tar.zst

# load allas tools
module load allas

# uncompress and remove the original files
zstdmt -d --rm $name"_R1_001.fastq.gz.zst
zstdmt -d --rm $name"_R2_001.fastq.gz.zst
```

### TO DO:

- DIAMOND blastx individually, concatenate all results
- METAXA pipeline (or outside snakemake?)
- unpack and remove original sequence files
- metadata file
- filter step for ResFinder mapping
- more log files as outputs eg. for ResFinder mapping
- more temp files to save memory
- add sample names to read headers in .bam files (eg. with samtools reheader)?
  or in fasta.gz files just before blastx?
- parse_diamond.py for CARD results
- metaphlan
- humann
-
