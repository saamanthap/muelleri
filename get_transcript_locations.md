Start with a file of transcript names, like this: 
```
TRINITY_DN51383_c0_g6_i2
TRINITY_DN51579_c0_g2_i4
TRINITY_DN52140_c7_g1_i12
TRINITY_DN51812_c3_g2_i10
TRINITY_DN57346_c1_g1_i3
```
Then use a for loop to grep the headers from your de novo transcriptome: 
```
for i in $(cat ./example_list_of_differentially_expressed_transcripts.txt); do grep -i "$i " ../muel/de_novo_assembly_trinity/muel_trinity_assembly_all_batches.Trinity.fasta >> names.txt;done
```
Now quickly delete the '>' signs. You can do this easily in vi:
```
vi names.txt
:%s/>//g
```
Now get the full fasta entries:
```
module load StdEnv/2020 seqtk/1.3
seqtk subseq ../muel/de_novo_assembly_trinity/muel_trinity_assembly_all_batches.Trinity.fasta names.txt > DE_genez.fasta
```
Now that you have the transcripts, you can blast them against any database you want. This is my script that performs a blast and then prints the hit with the highest bit score for each query: 
```
#!/bin/bash
#SBATCH --job-name=Blast
#SBATCH --cpus-per-task=12
#SBATCH --time=3:00:00
#SBATCH --mem=20GB
#SBATCH --output=Blast.%J.out
#SBATCH --error=Blast.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#run like this: sbatch blast.sh query_file blast_database output_name

module load StdEnv/2023 gcc/12.3 blast+/2.14.1

blastn -query ${1} -db ${2} -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -task dc-megablast |  sort -k1,1 -k12,12g -k13,13g | awk '!a[$1]++' > ${3}

#this script performs a blast search on a fasta file full of queries (usually these are taken from a DESeq2 analysis). The output is printed in a format that contains the information I want. The resulting file is sorted by the query name, then expect value, then bit score. For each queury sequence, the result with the highest bitscore is printed. This gives me a location for each differentially expressed transcript.
#before using this script, you'll have to make a query file using other commands (essentially, you take the transcript names that DESeq2 outputs and collect the fasta entries associated with them)
```
Now you want to join this with your DESeq2 results file so that you have the transcript location and expression info all on one line:
```
#!/bin/bash

#this script takes the filtered blast hits (meaning that each transcript has only one blast result), and adds information from deseq2 to the end of each line.
#it is really important that the blast result file is the FIRST argument, otherwise the column names will not make sense
# ${1} is the best blast results
# ${2} is the DESeq2 results
# ${3} is the name of the output file

touch ${3}
echo -e "transcript\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpajd" > ${3}

if [[ ! -f ${1} && ! -f ${2} ]]; then
        echo "Could not find input files. Run like this: ./join_blast_deseq2.sh blast_file deseq2_file out_file"
fi


join -t $'\t' <(sort -k 1 ${1}) <(sort -k 1 ${2}) >> ${3}
```
