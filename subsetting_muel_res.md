# Subsetting muel_res.txt
muel_res.txt contains the output of the DESeq2 analysis of differential expression. The only information that has been pre-filtered is transcripts with a count sum of less than 2 for all samples.  

The file looks a bit wonky, but that is because the tabs are different sizes. I have checked and everything is properly tab-delimited and can be reliable filtered by column: 
```transcriptID    baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
TRINITY_DN37947_c0_g1_i1        2.32589035223894        -1.10133580543906 1.0727632257983  -1.02663456292463       0.304592572718961       NA
TRINITY_DN816_c0_g1_i2  19.1745585951407        0.102368878209317       0.694342719096557  0.147432781238801       0.8827904325416 0.999986614612765
TRINITY_DN17376_c0_g1_i14       17.5829886832452        -1.60326581098384 1.74816592907749 -0.917113063649452      0.359083396197561       0.999986614612765
TRINITY_DN14109_c2_g2_i1        2.1227073860638 1.271017635168  1.28469026599311   0.989357255062146       0.322488380296208       NA
TRINITY_DN23737_c1_g4_i1        6.86482938282486        -0.1079001750607470.605669940780386        -0.178150124012628      0.858605078312968       0.999986614612765
```
Negative log2FoldChange values imply female-biased expression. Positive log2FoldChange values imply male-biased expression. I start by getting all the female-biased genes with an adjusted p-value < 0.05 and storing them in a new file:
```
awk -F'\t' '$3 < 0' muel_res.txt | awk -F'\t' '$7 < 0.05' > lfc_less_than_0_padj_less_than_0.05
```
Now I want to get all the transcript IDs and collect their sequences so that I can blast them against various databases. Here I prepare the query file:
```
cut -f 1 lfc_less_than_0_padj_less_than_0.05 | cat | grep -w -f - ../de_novo_assembly_trinity/muel_trinity_assembly_all_batches.Trinity.fasta | sed -e 's/>//g' > query.txt
```
Now use seqtk to extract the sequences for each entry: 
```
module load StdEnv/2023 seqtk/1.4
seqtk subseq ../de_novo_assembly_trinity/muel_trinity_assembly_all_batches.Trinity.fasta query.txt > fem_0.05.fasta
```
Now I can actually query a blast database. In this case, I want to get Xenbase gene IDs of the form "XBXL10_1g8966", so I'm going to query the Xenopus laevis reference genome. I had to make a database first:
```
makeblastdb -in XENLA_10.1_genome.fa -dbtype nucl -out XENLA_10.1_genome_blastable
```
Now I'm going to blast the fem_0.05.fasta file I prepared against the blast database I just made. The downside of using the X. laevis genome is that BLAST will find matches for exons, rather than whole transcripts or genes, since it is not splice-aware:
```
module load StdEnv/2023 gcc/12.3 blast+/2.14.1
blastn -query fem_0.05.fasta -db ../../2021_XL_v10_refgeno
me/XENLA_10.1_genome_blastable -outfmt 6 -out fem_0.05_to_XL_genome.out
```
