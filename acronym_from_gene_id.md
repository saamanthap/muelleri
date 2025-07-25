I want to make a file that contains xenbase GeneIDs/model IDs and gene acronyms, which can be used to understand my STAR DESeq2 results. I generated the STAR results using the XENLA_10.1_Xenbase.gtf file, which has entries that look like this: 
```
Chr2L   BestRefSeq      three_prime_utr 206122  206158  .
        -       .       gene_id "modelID:XBXL10_1g5996"; transcript_id "RefSeq:NM_001095632.1"; Dbxref "GeneID:733317" Dbxref "Genbank:NM_001095632.1"; ID "nbis-three_prime_utr-1307"; Note "The RefSeq transcript has 5 substitutions, 1 non-frameshifting indel and aligns at 99% coverage compared to this genomic sequence"; Parent "XBmRNA11210"; exception "annotated by transcript or proteomic data"; gbkey "mRNA"; gene "rts"; inference "similar to RNA sequence, mRNA (same species):RefSeq:NM_001095632.1"; original_biotype "three_prime_UTR"; partial "true"; previous_transcript_id "NM_001095632.1"; product "RECQL4-helicase-like protein";
Chr2L   Gnomon  gene    223300  234361  .       +       .           gene_id "Xenbase:XB-GENE-17340624"; Dbxref "GeneID:108707774" Dbxref "Xenbase:XB-GENE-17340624"; ID "XBXL10_1g5997"; Name "vps28.L"; Ontology_term "SO:0001217"; curie "Xenbase:XB-GENE-17340624"; gbkey "Gene"; gene "vps28.L"; gene_biotype "protein_coding";
```
Note that after "gene_id", some genes have a Xenbase gene ID, while others have a model ID. I want to collect an ID for each gene and match it to an acronym. In order to do this, I have to pull out the string inside *gene_id "  ";
```
