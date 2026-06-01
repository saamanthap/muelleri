## How to do a Trinity guided assembly
Instructions on the Trinity [Github page](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly). To perform a guided assembly, you need an input bamfile (the result of your RNASeq reads mapped to an appropriate reference genome). I'm going to align my reads again, since I have recently done a better job trimming them. I'll use STAR (which is splice-aware) and align to the borealis reference:
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
```
Then merge into a single file:
```
```
Now assemble:
```
```
