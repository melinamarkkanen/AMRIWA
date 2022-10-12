# Formulate output for rpoB counts into a table
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
