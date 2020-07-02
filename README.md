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

Batch job for post processing ResFinder mapping results in Puhti:
```
#!/bin/bash -l
#SBATCH --job-name=resfinder_results
#SBATCH --account=Project_2002265
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 2
#SBATCH --partition=small
#SBATCH --output=resfinder_results_out_%A_%a.txt
#SBATCH --error=resfinder_results_err_%A_%a.txt

module load biokit
cd /scratch/project_2002265/markkan5/AMRIWA/workflow/sorted_reads/

#name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2002265/markkan5/AMRIWA/sample_data.txt)                                ## somehow not correct

#echo -e $name > $name"_counts"                                                                                                 ## not working because the above is not correct
#samtools idxstats $name".bam" | grep -v "*" | cut -f3 >> ../resfinder_out/$name"_counts"                                       ## not working because the above is not correct

echo -e "GENE" > ../resfinder_out/gene_names
samtools idxstats BFH1_S123.bam | grep -v "*" | cut -f1 >> ../resfinder_out/gene_names

paste  ../resfinder_out/gene_names ../resfinder_out/*_counts > ../resfinder_out/ARG_genemat.txt
```
### TO DO:

- DIAMOND blastx individually, concatenate all results
- METAXA pipeline (or outside snakemake?)
- unpack and remove original sequence files
- metadata file
- filter step for ResFinder mapping
- more log files as outputs eg. for ResFinder mapping
- more temp files to save memory
- add sample names to read headers in .bam files? (eg. with samtools reheader)
- parse_diamond.py for CARD results
- metaphlan
- humann
-
