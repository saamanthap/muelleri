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
#!/bin/sh
#SBATCH --job-name=samtools_merge
#SBATCH --cpus-per-task=10
#SBATCH --time=8:00:00
#SBATCH --mem=50gb
#SBATCH --output=samtools_merge.%J.out
#SBATCH --error=samtools_merge.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load StdEnv/2023 samtools/1.22.1

bam_list=${1}
out_name=${2} #X_muelleri_alltads_merged.bam

samtools merge -b ${bam_list} -o ${out_name} --threads ${SLURM_CPUS_PER_TASK}
```
And sort:
```
#!/bin/sh
#SBATCH --job-name=samtools_sort
#SBATCH --time=4:00:00
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1
#SBATCH --output=samtools_sort.%J.out
#SBATCH --error=samtools_sort.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load StdEnv/2023 samtools/1.22.1

in=${1}

samtools sort -o ${in//.bam/_sorted.bam} ${in}

samtools index ${in//.bam/_sorted.bam}
```
Now assemble. You need to choose an appropriate value for 
```
#!/bin/sh
#SBATCH --job-name=trinity_guided_assembly
#SBATCH --cpus-per-task=32
#SBATCH --time=168:00:00
#SBATCH --mem=250G
#SBATCH --output=trinity_guided_assembly.%J.out
#SBATCH --error=trinity_guided_assembly.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2020
module load gcc/9.3.0 openmpi/4.0.3
module load trinity/2.14.0 samtools/1.17 jellyfish/2.3.0
module load salmon/1.4.0
module load python/3
module load scipy-stack/2023a

bam=${1}

Trinity --genome_guided_bam ${bam} \
        --genome_guided_max_intron 20000 \
        --max_memory 200G \
        --CPU ${SLURM_CPUS_PER_TASK} \
        --output /home/samp/projects/rrg-ben/for_Sam/muel/sam_trinity_assembly/trinity_guided_assembly/
```
