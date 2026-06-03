## How to do a Trinity guided assembly
Instructions on the Trinity [Github page](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly). To perform a guided assembly, you need an input bamfile (the result of your RNASeq reads mapped to an appropriate reference genome). Guided transcriptome assembly is only a good idea if you have a relatively complete reference genome (a highly fragmented reference genome will not work.) I'm going to align my reads again, since I have recently done a better job trimming them. I'll use STAR (which is splice-aware) and align to the borealis reference:
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --array=0-9
#SBATCH --cpus-per-task=10
#SBATCH --time=6:00:00
#SBATCH --mem=100gb
#SBATCH --output=STAR_map.%J.%a.out
#SBATCH --error=STAR_map.%J.%a.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -euxo pipefail

module load StdEnv/2023 star/2.7.11b

# run like this: sbatch wip_STAR
# run from the directory you want output and error files to save in

dir=/home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed/
declare -a forward=(${dir}/*R1_bbduk.fastq.gz)
fread=${forward[$SLURM_ARRAY_TASK_ID]}
rread=${fread//R1_bbduk.fastq.gz/R2_bbduk.fastq.gz}
outdir=/home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed_star_mapped/
outpre=$(basename $fread R2_bbduk_fastqc.zip)
ref=/home/samp/projects/rrg-ben/for_Sam/Austin_genome_Xborealis/
gtf=/home/samp/projects/rrg-ben/for_Sam/Austin_genome_Xborealis/Xbo.cds.gff

STAR --genomeDir ${ref} \
--sjdbGTFfile ${gtf} \
--runThreadN ${SLURM_CPUS_PER_TASK} \
--readFilesIn ${fread} ${rread} \
--outFileNamePrefix ${outdir}${outpre} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--alignEndsProtrude 10 ConcordantPair \
--outFilterMatchNminOverLread 0.5 \
        --readFilesCommand zcat
```
Since Trinity only accepts a single bam file as input, you need to add readgroups, then merge the bams into one before you can actually do the assembly. Start by adding ReadGroups:
```
#!/bin/bash
#SBATCH --job-name=readgroups
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=24gb
#SBATCH --output=readgroups.%J.%a.out
#SBATCH --error=readgroups.%J.%a.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2023 gatk/4.6.1.0


#load bams into array (run from the folder than contains bams)
bams=( ./*Aligned.sortedByCoord.out.bam )


in=${bams[$SLURM_ARRAY_TASK_ID]}
sample_name=$(basename $input_file Aligned.sortedByCoord.out.bam)

module load StdEnv/2023  picard/3.1.0

gatk AddOrReplaceReadGroups \
        -I ${in} \
        -O ${sample_name}_rg.bam \
        -RGID ID_${sample_name} \
        -RGLB lib_${sample_name} \
        -RGPL "Illumina" \
        -RGPU unit_${sample_name} \
        -RGSM ${sample_name} \
        --CREATE_INDEX true
```
Then merge into a single file:
```
```
Now assemble:
```
```
