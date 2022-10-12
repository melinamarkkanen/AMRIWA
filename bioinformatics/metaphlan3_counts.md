# Metaphlan3 to obtain counts
```
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./sample_names.txt)

metaphlan $name"_R1_trimmed.fastq.gz" \
        --nproc $SLURM_CPUS_PER_TASK \
        --bowtie2out $name".bowtie2.bz2" \
        --input_type fastq > $name"_profile.txt" \
        --bowtie2db metaphlan_databases/ \
        -t rel_ab_w_read_stats
```

# Gathering and wrangling the output and the tax table
```
./merge_metaphlan_tables_abs.py *_profile.txt > metaphlan3_counts_merged_abundance_table.txt
sed -i 's/_profile//g' metaphlan3_counts_merged_abundance_table.txt

sed -n '2p' metaphlan3_counts_merged_abundance_table.txt > merged_counts_abundance_table_genus.txt
grep -E "g__" metaphlan3_counts_merged_abundance_table.txt >> merged_counts_abundance_table_genus.txt
tr '|' ';' <merged_counts_abundance_table_genus.txt > merged_metaphlan_counts_genus.txt

# tax table
awk '{print $1}' merged_metaphlan_counts_genus.txt > tax_table_count_metaphlan3_genus
grep -E "g__" merged_metaphlan_counts_genus.txt > tax_table_count_metaphlan3_genus.txt
```


#### (Run similarly for the relative abundances with the `merge_metaphlan_tables.py` program)
