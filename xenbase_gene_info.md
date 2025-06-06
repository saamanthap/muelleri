# Getting gene info from 
First download databases of gene info from Xenbase. I chose a file of general info about genes, gene expression data in laevis, and stage amd tissue expression profiles:
```
wget https://download.xenbase.org/xenbase/GenePageReports/GenePageGeneralInfo_ManuallyCurated.txt
wget https://download.xenbase.org/xenbase/GenePageReports/GeneExpression_laevis.txt
wget https://download.xenbase.org/xenbase/GenePageReports/GenePageAnatomyOntologyMapping.txt
```
Have a .txt file of gene names formatted with one gene per line. These were pulled from laevis transcriptome blast results:
```
LOC108715912
LOC108718162
LOC108718162
aars1.L
bchel.L
copb1.L
copb1.L
```
Use the file of gene names as an input for grep: 
NEED TO FIX THIS SO THAT IT ONLY MATCHES WHOLE WORDS
```
grep -f ../XENLA_10.1_Xenbase.transcripts/gene_names.txt GeneExpression_laevis.txt
```
GenePageAnatomyOntologyMapping.txt and GenPageGeneralInfo_ManuallyCurated.txt do not contain .S or .L gene names, so you have to alter the gene names:
NEED TO FIX THIS SO THAT IT ONLY MATCHES WHOLE WORDS
```
cut -d'.' -f1 ../XENLA_10.1_Xenbase.transcripts/gene_names.txt | grep -F - ./GenePageAnatomyOntologyMapping.txt | less
```
If you need to eliminate identical gene names from the list:  
*Note: -f- allows grep to take standard input*
```
 awk '{a[$1]++} END{for(b in a) print b}' ../XENLA_10.1_Xenbase.transcripts/gene
_names.txt | grep -f- ./GeneExpression_laevis.txt | less
```
