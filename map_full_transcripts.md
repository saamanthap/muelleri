I want to use STAR to map full-length transcripts (assembled from X. muelleri RNA-Seq reads) to the *X. borealis* reference genome. Previously, I tried to tackle this problem using pblat (see other code in this repo), but it takes SO LONG, since I have 500000 transcripts in the assembly.
```
STAR --genomeDir /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index/ \
--runThreadN 6 \
--readFilesIn Mov10_oe_1.subset.fq \
--outFileNamePrefix ../results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
```
