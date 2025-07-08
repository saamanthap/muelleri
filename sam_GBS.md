# Sam_GBS.md
### Modified versions of Ben's scripts from [2020_GBS](https://github.com/evansbenj/2020_GBS.md) which can be run in parallel on Brian's machine.

## Trimmomatic
Trimmomatic is used to trim adapters and primers that remain on reads from library prep and sequencing. Make sure that you have an adaptor file in fasta format that contains the sequences of all the adaptors and primers you want to trim.  

My adaptor file is called **TruSeq2_and_3-PE-2.fa**. It looks like this: 
```
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```
My script is called **trim.Sam** and I store it on the scratch drives on Brian's machine:
```
#!/bin/sh

# run like this: parallel -j x ../trim.Sam ::: *_R1_001.fastq.gz &
# where x is the maximum number of jobs to run at once
# make sure to run from the directory that contains the untrimmed reads

adapt=/2/scratch/TruSeq2_and_3-PE-2.fa  #full path to the adaptor file
trim=/opt/local/trimmomatic/trimmomatic-0.39.jar  #full path to the jar file on Brian's machine
base=$(basename $1 _R1_001.fastq.gz) #removing the suffix from the file name
forward=$1 #input forward read
reverse=${base}_R2_001.fastq.gz #input reverse read

#making sure the input files exist
if [ ! -f $forward ]; then
        echo "Error: R1 not found" 2> ${base}_err.txt
        exit 1
fi

if [ ! -f $reverse ]; then
        echo "Error: R2 not found" 2>> ${base}_err.txt
        exit 1
fi

#trimming the reads
java -jar $trim PE \
        -trimlog ${base}_log.txt \
        $forward \
        $reverse \
        ${base}_trim_R1.fq.gz \
        ${base}_trim_single_R1.fq.gz \
        ${base}_trim_R2.fq.gz \
        ${base}_trim_single_R2.fq.gz \
        ILLUMINACLIP:$adapt:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
## FastQC
Use fastqc to generate reports on the quality of your reads. It is prudent to generate reports before and after trimming to compare.   

For a single file: 
```
fastqc X_muelleri_tad31_S11_L001__trim_R1.fq.gz*
```
To quickly generate many reports at once: (x is the maximum number of reports to generate at once... I have learned that 20 reports at once consumes all the swap memory on info114)
```
parallel -j x fastqc ::: *trim_R?.fq.gz
```
For the above command I used a regular expression that will catch all the files I want to generate fastqc reports for and feed them as input to parallel.  

Sometimes I have trouble on Brian's machine because commands like fastqc are not found. If this happens, you can either a) temporarily add fastqc to your path:
```
export PATH="/path/to/directory/containing/executable:$PATH"
```
Or b) navigate to home, and add the above line to the end of your .bash_profile file.
## Preparing the reference genome
First you need a reference genome! You can download the X. laevis ref genome from Xenbase. Get the .fa file and the .gtf file:
```
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/latest/XENLA_10.1_Xenbase.gtf.gz
```
You can unzip the files if you want (though you don't need to):
```
gunzip file_name
```
Next you'll need to index the genome:
```
bwa index /2/scratch/samp/XENLA_10.1_genome/XENLA_10.1_genome.fa
```
The indexing output is several files: .amb .ann .bgz .bwt .pac .sa  

Next generate a .fai file:
```
samtools faidx XENLA_10.1_genome.fa
```
Then a .dict file. Since you're using picard, which runs on java, you call the program a little bit differently:
```
java -jar /opt/local/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=XENLA_10.1_genome.fa OUTPUT=XENLA_10.1_genome.dict
```
## Map data
Next you'll use the bwa mem command to map your reads to the reference genome. This is a good process to parallelize. My script is called **for.Sam** (thanks to Brian for his help!):
```
#!/bin/bash

ref='/2/scratch/samp/XENLA_10.1_genome/XENLA_10.1_genome.fa'
# NOTE the ref sequence has to be indexed FIRST

rf=$1
rr=`basename $1 _R1.fq.gz`_R2.fq.gz
sortedbam=`basename $rf _R1.fq.gz`_sorted.bam
# echo "bwa mem $ref $rf $rr -t 5 | samtools view -Shu - | samtools sort - -o $sortedbam"
# echo "samtools index $sortedbam"

bwa mem $ref $rf $rr -t 5 | samtools view -Shu - | samtools sort - -o $sortedbam
samtools index $sortedbam

# Run as: parallel -j 5 for.Sam *_R1.fq.gz &>> bwa_err.txt &
# limit to five (or whatever) jobs simultaneously
# Run from the directory containing trimmed fastq files and give the correct relative path to for.Sam
```
## Add readgroups
My script is called **readgroups.Sam**:
```
#!/bin/sh

# Run using parallel -j x readgroups.Sam ::: *_sorted.bam &>> rg_err.txt &
# Run from the directory that contains the sorted.bam files
# Make sure that you give the correct relative (or absolute path) to the readgroups.Sam file

picard=/opt/local/picard-tools/picard.jar
input=$1
output=$(basename $1 _sorted.bam)_rg.bam
rglibrary=$(basename $1)
rgsample=$(basename $1)
rgunit=$(basename $1)

java -jar $picard AddOrReplaceReadGroups \
        I=$input \
        O=$output \
        RGID=4 \
        RGLB=$rglibrary \
        RGPL=ILLUMINA \
        RGPU=$rgunit \
        RGSM=$rgsample
```
## GATK HaplotypeCaller
My script is called **haplo.Sam**:
```
#!/bin/sh

# Run like this: parallel -j x ../haplo.Sam ::: *_rg.bam &>> haplo_err.txt &
# x is the maximum number of jobs to run at one time
# Run from the directory that contains the _rg.bam files - in my case, \2\scratch\samp\muel_bam\
# Make sure to give the correct relative path to haplo.Sam
# In this case I am using the laevis reference genome (downloaded from xenbase)

gatk=/home/samp/opt/local/gatk/GenomeAnalysisTK.jar
ref=/2/scratch/samp/XENLA_10.1_genome/XENLA_10.1_genome.fa
input=$1
output=$(basename $1 _rg.bam).g.vcf.gz

java -Xmx24G -jar $gatk HaplotypeCaller \
        -R $ref \
        -I $input \
        -O $ouput \
        -ERC GVCF
```

