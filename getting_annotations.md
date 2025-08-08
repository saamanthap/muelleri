Get the acronym after GenePageIDs: 
```
awk '{print gensub(/.*GenePageIDs:([^|]+)\|.*/, "\\1", "g")}' m_f_laevis_transcriptome_blast_results.out > m_f_annotations
```
Find a string inside () and replace the entire line with the string. This is a command for inside vi: 
```
%s/.*(\(.*\)).*/\1/
```
I want to add LFC and padj to the end of my modified blast file: (NEED TO EXPLAIN HOW + WHY I MADE THIS BLAST FILE)
```
awk -v OFS="\t" 'NR==FNR { a[$1]=$3" "$7; next } $1 in a { print $0, a[$1] }' deseq2_results blast_results
```
Explanation: 
- Make an associative array using the first column of the deseq2_results file as the key (transcript names). Inside the array I'm going to store the LFC and the padj, which are column $3 and $7 respectively. This way I can use the transcript 
