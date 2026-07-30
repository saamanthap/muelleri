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
Now use BLAST to get corresponding regions in borealis. I set the max target seqs to 5 arbitrarily because I'm going to end up filtering for the best hit anyway.
```
blastn -task dc-megablast -query ../2021_XL_v10_refgenome/XENLA_10.1_GCF_exononly.fasta -db Xbo.v1.fa_blastable -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" -out laevis_exons_to_borealis.out -max_target_seqs 5
```
