# AMRIWA
## AMRIWA metagenome analysis

### Analysis using the Snakemake workflow

#### Create virtual environment
```
conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```

#### The Snakemake file and the environment `.yml` files are in the `workflow` folder.  

#### Dry-run
```
snakemake --use-conda -np
```
#### Run in batch job
```
snakemake --use-conda -j $SLURM_CPUS_PER_TASK
```

# Data wrangling

## Tidy Metaphlan3 outputs prior loading into R
```
sed -i 's/_profile//g' merged_abundance_table.txt
sed -n '2p'  merged_abundance_table.txt > merged_abundance_table_species.txt
sed -i 's/_profile//g' merged_abundance_table_species.txt
grep -E "s__" merged_abundance_table.txt >> merged_abundance_table_species.txt
tr '|' ';' <merged_abundance_table_species.txt > mod_merged_abundance_table_species.txt

# Tax table
awk '{print $1}' mod_merged_abundance_table_species.txt > tax_table_metaphlan
sed '1d' -i tax_table_metaphlan
```

## Tidy MGE mapping outputs prior loading into R
```
# OTU table
cp MGE_genemat.txt cp_MGE_genemat.txt
sed -i -e 's/\[//g' cp_MGE_genemat.txt
sed -i -e 's/\]//g' cp_MGE_genemat.txt
sed -i -e 's/\"//g' cp_MGE_genemat.txt

# Tax table
# Three accessions are incorrectly formatted due to the quotation marks, 
# those have to bee modified manually.
sed -i -e 's/\[//g' cp_MGE.fasta
sed -i -e 's/\]//g' cp_MGE.fasta
sed -i -e 's/\"//g' cp_MGE.fasta
```
## Create the tax table for ARGs with ARG clusters using CD-HIT (v4.8.1) with 90 % identity
```
cd-hit-est -i resfinder.fasta -o clusters_resfinder.fasta -c 0.90

cat resfinder_db/phenotypes.txt | awk '{ print $1, $2}' > genes_classes.txt
sed -i '1d' genes_classes.txt
sed -i 's/ /\t/g' genes_classes.txt

# Name all variants in the same cluster with a common name (the one with *) -> cluster_names.txt
sed -i 's/;/\t/g' cluster_names.txt
# Get all ARG names
sed -i 's/>//g' ARG_genes.txt
grep -f ARG_genes.txt cluster_names.txt > clusters_tax_table.txt
```

## Formulate output for rpoB counts into a table. Include counts for both reads here.
### However, only forward (R1) reads will be used in the analysis (See R scripts).
```
ls *_hmm_out.txt | sed 's/_hmm_out.txt//g' > HMM_names.txt
cat *_HMM_count.txt > HMM_counts.txt

echo -e "sample_name" > HMM_sample_names.txt
awk '0 == (NR + 1) % 2'  HMM_names.txt | sed 's/_R1//g' >> HMM_sample_names.txt

# read1
echo -e "R1" > R1_HMM_counts
awk '0 == (NR + 1) % 2'  HMM_counts.txt >> R1_HMM_counts

# read2
echo -e "R2" > R2_HMM_counts
awk '0 == NR % 2'  HMM_counts.txt >> R2_HMM_counts

paste -d"\t" HMM_sample_names.txt R1_HMM_counts R2_HMM_counts > HMM_RESULT_TABLE.txt
```

# Metaxa2 (v.2.2.1)
## Run in batch job
```
metaxa2 -1 ~/AMRIWA/workflow/trimmed_data/$name"_R1_trimmed.fastq.gz" \
        -2 ~/AMRIWA/workflow/trimmed_data/$name"_R2_trimmed.fastq.gz" \
        -f fastq -z gzip -o $name --align none --graphical F --cpu $SLURM_CPUS_PER_TASK --plus

metaxa2_ttt -i $name".taxonomy.txt" -o $name

# Combine outputs to get genus level taxa
metaxa2_dc -o metaxa_genus.txt *level_6.txt
```

# Get gene lengths for ARGs and MGEs using SeqKit (v0.12.1)
```
seqkit fx2tab --length --name --header-line resfinder.fasta > resfinder_lengths.txt
seqkit fx2tab --length --name --header-line MGE.fasta > MGE_lengths.txt
```

# Assembly for a subset of metagenomic samples using MEGAHIT (v1.2.8)
and MetaQUAST (QUAST v5.0.2) for the quality control. 
## Run in batch job.
### (BH02, BH10, BH13, BH48, BH52, BFH19, BFH26, BFH27, BFH29, BFH35, BFH42, FH1, FH7, FH9)
```
megahit -1 $name"_R1_trimmed.fastq.gz"  -2 $name"_R2_trimmed.fastq.gz"  \
         -o $name"_assembly" -t $SLURM_CPUS_PER_TASK --min-contig-len 1000 -m 32000000000

# MetaQUAST
cd $name"_assembly"
metaquast.py -t $SLURM_CPUS_PER_TASK --no-plots -o $name"_assembly_QC" final.contigs.fa
```

# Run BLAST search between the obtained contigs and ResFinder database
# in order to identify the contigs encoding plasmid mediated colistin ARGs (mcr genes).
```
cd $name"_assembly"
blastn -query final.contigs.fa \
        -subject ~/AMRIWA/workflow/resfinder_db/resfinder.fasta \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp slen qstart qend qlen" \
        -out resfinder_blast_out.txt -max_target_seqs 50 -perc_identity 90
``` 

## Extract contigs encoding the ARG of interest (mcr gene) using SeqKit tools (v0.12.1).
```
seqkit grep -p <contigID> final.contigs.fa > <sample>_CONTIG_<contigID>.fasta
```

## Visualise assembled contigs using Bandage (https://github.com/rrwick/Bandage.git)
### Prepare fastg files for visualisation.
```
megahit_toolkit contig2fastg 141 k141.contigs.fa > $name".fastg"
```

### Create BLAST database for Bandage
```
grep "in" MGE.fasta > int_qac_ID_list.txt
grep "qacE" MGE.fasta >> int_qac_ID_list.txt

sed -i 's/>//' int_qac_ID_list.txt
seqkit grep -n -f int_qac_ID_list.txt MGE.fasta > int_qac.fasta

cat resfinder.fasta int_qac.fasta > Bandage_BLAST_db
```

### Run locally.
```
Bandage.app/Contents/MacOS/Bandage load $name".fastg" --draw --query Bandage_BLAST_db
```
