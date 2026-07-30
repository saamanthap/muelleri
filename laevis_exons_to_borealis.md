Pull out exon entries from gff
```
awk -F'\t' '$3 ~ /exon/' XENLA_10.1_GCF.gff3 > XENLA_10.1_GCF_exononly.gff3
```
Reformat lines
```
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $4, $5, $9, ".", $7}' TEMP_XENLA_10.1_GCF_exononly.gff3 > XENLA_10.1_GCF_exononly.gff3
```
Use bedtools getfasta to extract sequences
```
bedtools getfasta -fi XENLA_10.1_genome.fa -bed XENLA_10.1_GCF_exononly.gff3 -name -fo XENLA_10.1_GCF_exononly.fasta
```
Now use BLAST to get corresponding regions in boralis
```

```
