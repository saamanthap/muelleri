Make the samples.txt file: 
```
for r1_file in *_R1.fq.gz; do base_name="${r1_file%_R1.fq.gz}"; r2_file="${base_name}_R2.fq.gz"; if [[ -f "$r2_file" ]]; then common_prefix=$(echo "$base_name" | cut -d'_' -f1-3); species_prefix=$(echo "$base_name" | cut -d'_' -f1-2); echo -e "${species_prefix}\t${common_prefix}\t${r1_file}\t${r2_file}"; fi; done > processed_filenames.txt
```
The text file looks like this: 
```
X_muelleri      X_muelleri_tad31        X_muelleri_tad31_S11_L001__trim_cut_polyA_R1.fq.gz      X_muelleri_tad31_S11_L001__trim_cut_polyA_R2.fq.gz
X_muelleri      X_muelleri_tad32        X_muelleri_tad32_S12_L001__trim_cut_polyA_R1.fq.gz      X_muelleri_tad32_S12_L001__trim_cut_polyA_R2.fq.gz
```
