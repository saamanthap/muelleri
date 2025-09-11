# Sam_GBS.md
### Modified versions of Ben's scripts from [2020_GBS](https://github.com/evansbenj/2020_GBS.md), some of which have been modified for Brian's machines, and others that work on Nibi.

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
Version of bwa mem script made for nibi (called **wip_bwamem**)
```
#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --cpus-per-task=4
#SBATCH --array 0-9
#SBATCH --time=10:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.%j.out
#SBATCH --error=bwa_align.%J.%j.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#this script is used to map trimmed reads to a reference genome. the output SAM file is piped to samtools view, which reads the file. finally, the read SAM file is piped to samtools sort, which outputs the final sorted BAM file.

#run like this: sbatch wip_bwamem path_to_ref path_to_fq_files

module load StdEnv/2023 bwa/0.7.18 gcc/12.3 samtools/1.20

ref=${1}
fq_dir=${2}

declare -a fwd=(${fq_dir}*trim_cut_polyA_R1.fq.gz)

current_fwd=${fwd[${SLURM_ARRAY_TASK_ID}]}
current_rev=${current_fwd/R1.fq.gz/R2.fq.gz}
sorted=${current_fwd/trim_cut_polyA_R1.fq.gz/sorted.bam}

if [[ -e ${current_fwd} && -e ${current_rev} ]]; then
        bwa mem ${ref} ${current_fwd} ${current_rev} -t ${SLURM_CPUS_PER_TASK} | samtools view -Shu - | samtools sort - -o ${sorted}
        samtools index ${sorted}
fi


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
I also have a version of haplotypeCaller for use on nibi (I reccomend using this instead of working on Brian's machine, on which big jobs occaisionally fail). This script is called **
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --array=0
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH --output=HaplotypeCaller.%A.%a.out
#SBATCH --error=HaplotypeCaller.%A.%a.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# This script will read in the *trim_rg.bam file names in a directory, and make and execute the GATK command "HaplotypeCaller" on these files.
# Run like this: sbatch path/to/wip_haplo *file_ending

module load StdEnv/2023 gatk/4.6.1.0

ref=/home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/XENLA_10.1_genome.fa
in=/home/samp/projects/rrg-ben/for_Sam/muel/muel_bam/${1}
outdir=/home/samp/projects/rrg-ben/for_Sam/muel/new_muel_sorted

declare -a file=(${in})

current=${file[${SLURM_ARRAY_TASK_ID}]}
base=$(basename ${current} .bam)

mkdir -p ${outdir}

gatk --java-options -Xmx24G HaplotypeCaller -I ${current} -R ${ref} -O ${outdir}/${base}.g.vcf -ERC GVCF --native-pair-hmm-threads 4

```
## Combine GVCFs 
My script is called **wip_combineGVCF**. You'll have to pick an interval to combine the GVCFs across (in this case I chose Chr4L). If you want to combine GVCFs for all chromosomes, use the array version of this script, which is below.
```
#!/bin/sh
#SBATCH --job-name=GenomicsDBImport
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --time=60:00:00
#SBATCH --mem=100gb
#SBATCH --output=GenomicsDBImport.%J.out
#SBATCH --error=GenomicsDBImport.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# This script will read in the *.g.vcf file names in a directory, and
# make and execute the GATK command "GenotypeGVCFs" on these files.

# execute like this:
# sbatch 2021_GenomicsDBImport.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/vcf/ chr1L temp_path_and_dir db_path_and_dir_prefix

# argument 1 is the ref genome
# argument 2 is the .g.vcf files
# argument 3 is the chromosome
# argument 4 is the path and dir for temp files
# argument 5 is the path to the database and the chromosome

# (sam's changes) run like this: sbatch wip_combineGVCF
# this command outputs a GenomicsDB workspace

# pretty sure i don't need to load nixpkgs, since Graham already has gatk
module load StdEnv/2023 gatk/4.6.1.0

ref=/home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/XENLA_10.1_genome.fa
vcfdir=/home/samp/projects/rrg-ben/for_Sam/muel/new_muel_sorted/
chrom=Chr4L
workspace=/home/samp/projects/rrg-ben/for_Sam/muel/DBI_combineGVCF
temp=/home/samp/projects/rrg-ben/for_Sam/muel/DBI_combineGVCF_temp/


commandline="gatk --java-options -Xmx20G GenomicsDBImport -R $ref"
for file in ${vcfdir}*rg.g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${chrom} --tmp-dir $temp --batch-size 50 --genomicsdb-workspace-path ${workspace}/ --max-num-intervals-to-import-in-parallel ${SLURM_CPUS_PER_TASK} --genomicsdb-shared-posixfs-optimizations"

${commandline}

```
CombineGVCFs for all chromosomes (submit each chromosome separately, in an array). Here, I hard-coded in the names of all the chromosomes. This script is called **wip_combineGVCF_array**
```
#!/bin/sh
#SBATCH --job-name=GenomicsDBImport
#SBATCH --array=0-17
#SBATCH --cpus-per-task=10
#SBATCH --time=6:00:00
#SBATCH --mem=24gb
#SBATCH --output=GenomicsDBImport.%J.%j.out
#SBATCH --error=GenomicsDBImport.%J.%j.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# This script will read in the *.g.vcf file names in a directory, and
# make and execute the GATK command "GenotypeGVCFs" on these files.

# (sam's changes) run like this: sbatch wip_combineGVCF
# this command outputs a GenomicsDB workspace

# pretty sure i don't need to load nixpkgs, since Graham already has gatk
module load StdEnv/2023 gatk/4.6.1.0

# i'm making this script into an array to run each chromosome as a separate job

declare -a chroms=("Chr1L" "Chr1S" "Chr2L" "Chr2S" "Chr3L" "Chr3S" "Chr4L" "Chr4S" "Chr5L" "Chr5S" "Chr6L" "Chr6S" "Chr7L" "Chr7S" "Chr8L" "Chr8S" "Chr9_10L" "Chr9_10S")

ref=/home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/XENLA_10.1_genome.fa
vcfdir=/home/samp/projects/rrg-ben/for_Sam/muel/new_muel_sorted/
chrom=${chroms[${SLURM_ARRAY_TASK_ID}]} # add coordinates of sex-specific region
workspace=/home/samp/projects/rrg-ben/for_Sam/muel/DBI_combineGVCF_all_chroms
temp=/home/samp/projects/rrg-ben/for_Sam/muel/DBI_combineGVCF_all_chroms/DBI_combineGVCF_temp/

commandline="gatk --java-options -Xmx20G GenomicsDBImport -R $ref"
for file in ${vcfdir}*rg.g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${chrom} --tmp-dir $temp --batch-size 50 --genomicsdb-workspace-path ${workspace}_${chrom}/ --reader-threads ${SLURM_CPUS_PER_TASK}"

${commandline}

```
## GenotypeGVCFs
Since I used GenomicsDBImport for the last step, I must use the -V gendb:// flag. Also note that you MUST use the same version of gatk for this and the previous step.
```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=23gb
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# This script will read in the *.g.vcf file names in a directory, and
# make and execute the GATK command "GenotypeGVCFs" on these files.

# execute like this:
# sbatch 2021_GenotypeGVCFs_DB.sh

ref=/home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/XENLA_10.1_genome.fa
# in this case, db will be ../DBI_combineGVCF_Chr4L (if i run from the genotype_gvcf directory)

module load StdEnv/2023 gatk/4.6.1.0

workspace=/project/6019307/for_Sam/muel/DBI_combineGVCF_Chr4L
db=$(basename ${workspace})

commandline="gatk --java-options -Xmx18G GenotypeGVCFs -R ${ref} -V gendb://${workspace} -O ${db}_out.vcf"

${commandline}
```
I skipped VariantFiltration and SelectVariants... this means that my data is noisier. In the future, it is worth using these to clean up less reliable variants. (See Ben's GBS page for more info.)  
## vcftools
I used vcftools to output a new vcf file that contains SNPs only. You can run this from the head, and should only take a few seconds:
```
vcftools --vcf input_file.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only
```
After these steps are done, use my [python script](https://github.com/saamanthap/muelleri/blob/main/python_find_heterozygosity.md) for flagging heterozygous or sex-specific sites. 


