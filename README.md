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
Download  data from allas:

First set password for Allas for project. Must be run every time prior to actual batch job.
```
module load allas
allas-conf -k project_2002265
```

Then download data with batch job:
```
#!/bin/bash
#SBATCH -J download_data
#SBATCH --account=project_2002265
#SBATCH -o download_data_%j_out.txt
#SBATCH -e download_data_%j_err.txt
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -p small

# go to target dir
cd /scratch/project_2002265/markkan5/AMRIWA/workflow/data

# make sure the connection to Allas is open
source /appl/opt/allas-cli-utils/allas_conf -f -k $OS_PROJECT_NAME

# download the whole bucket
swift download 2002265_Melina_Markkanen_AMRIWA_metagenomes

# remove the files that are not used here
rm *_ameta
rm FASTQC.tar.zst
```
Decompress data:
```
#!/bin/bash
#SBATCH -J decompress
#SBATCH --account=project_2002265
#SBATCH -o decompress_%j_out.txt
#SBATCH -e decompress_%j_err.txt
#SBATCH -t 10:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -p small

# go to dir and load zstd tool
cd /scratch/project_2002265/markkan5/AMRIWA/workflow/data
module load allas

# decompress reads R1 and R2
zstd -d --rm -q *_R1_001.fastq.gz.zst
zstd -d --rm -q *_R2_001.fastq.gz.zst
```

Install Metaxa2:
```
wget https://microbiology.se/sw/Metaxa2_2.2.1.tar.gz
tar -zxvf Metaxa2_2.2.1.tar.gz
./install_metaxa2
```

Running Metaxa2 in Puhti:
```
#!/bin/bash -l
#SBATCH -J metaxa
#SBATCH --account=project_2002265
#SBATCH -o metaxa2_out_%A_%a.txt
#SBATCH -e metaxa2_err_%A_%a.txt
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-96
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH -p small

cd /scratch/project_2002265/markkan5/AMRIWA/Metaxa2

module load biokit

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)

metaxa2 -1 /scratch/project_2002265/markkan5/AMRIWA/workflow/trimmed_data/$name"_R1_trimmed.fastq.gz" \
        -2 /scratch/project_2002265/markkan5/AMRIWA/workflow/trimmed_data/$name"_R2_trimmed.fastq.gz" \
        -f fastq -z gzip -o $name --align none --graphical F --cpu $SLURM_CPUS_PER_TASK --plus

metaxa2_ttt -i $name".taxonomy.txt" -o $name
```

### TO DO:

- metaxa2 pipeline -> or outside snakemake
-  
`
