# Muelleri Project Diary
## Counts  
* [Automatically collect abundance files from a directory to compile counts](https://github.com/saamanthap/muelleri/blob/main/sam_2025_compile_kallisto_muel.md)

## DE analysis in R  
[DESeq2 workflow](https://github.com/saamanthap/muelleri/blob/main/muelleri_differential_gene_expression.md)
* Filtering low counts, normalization, volcano plot, PCA, parallel coordinates plot   
* Collect transcripts with male or female biased expression (negative or positive LFC)   
* Collect transcripts with expression in females only  
* Collect transcripts with expression in females only with LFC <= -20


## DE analysis on command line  
[Useful bash commands](https://github.com/saamanthap/Ideas_for_Jade.md/blob/main/2025_checking_DE_transcript_locations.md)  
#### For female-only transcripts with LFC <= -20... strategies to identify genes of interest: 
* Get a list of transcripts that blast to Chr4L 110000000-147000000 on the X. laevis genome
* Get a list of transcripts that match the X. laevis transcriptome with alignment length/query length >= 0.85 (in order to identify "true" matches to protein-coding transcripts)
* Pull out unique transcript names based on the outfile format (useful if a transcript finds multiple matches, but can be misleading for matches with low alen/qlen ratio)
#### Strategies to get annotations: 
*Note: XENLA_10.1_Xenbase.transcripts has conveniently formatted gene names but 2025_XL_transcriptome does not*
* Collect fasta files from a list of names of transcripts of interest. Blast against the human genome (or another appropriate annotated genome), making sure to use the -task blastn flag
* If you have gene IDs and need gene names, you can search a line-separated list on the [Gene Ontology](https://www.pantherdb.org/) website


## For the future
#### Transcripts of interest
* Excluding transcripts with zero expression in males excludes transcripts that might have a female-specific insertion, and therefore still have *some* expression in males --> try looking at female-biased transcripts with a low level of expression in males
* Sex could be determined via haploinsuffciency, where males have two functional alleles, and females have only one. This kind of transcript would have roughly twice as much expression in males compared to females --> try looking at transcripts with male-biased expression, especially those located in sex-specific region on Chr4L

