```
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./sample_names.txt)

metaphlan $name"_R1_trimmed.fastq.gz" \
        --nproc $SLURM_CPUS_PER_TASK \
        --bowtie2out $name".bowtie2.bz2" \
        --input_type fastq > $name"_profile.txt" \
        --bowtie2db metaphlan_databases/ \
        -t rel_ab_w_read_stats
```
