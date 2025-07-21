Make the samples.txt file: 
```
for r1_file in *_R1.fq.gz; do base_name="${r1_file%_R1.fq.gz}"; r2_file="${base_name}_R2.fq.gz"; if [[ -f "$r2_file" ]]; then common_prefix=$(echo "$base_name" | cut -d'_' -f1-3); echo -e "${common_prefix}\t${r1_file}\t${r2_file}"; fi; done > processed_filenames.txt
```
