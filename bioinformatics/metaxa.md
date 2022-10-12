# metaxa2
```
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../sample_names.txt)

metaxa2 -1 $name"_R1_trimmed.fastq.gz" \
        -2 $name"_R2_trimmed.fastq.gz" \
        -f fastq -z gzip -o $name --align none --graphical F --cpu $SLURM_CPUS_PER_TASK --plus

metaxa2_ttt -i $name".taxonomy.txt" -o $name
```

# metaxa_dc
```
metaxa2_dc -o metaxa_genus.txt *level_6.txt
```
