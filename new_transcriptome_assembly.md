# New and improved transcriptome assembly!
I am trying to assemble a new transcriptome where I have less independently assembled alleles. This will increase my statistical power when it comes to the differential expression analysis. The idea here is that I will trim more thoroughly (using bbduk instead of trimmomatic and cutadapt), and then assemble a male and female transcriptome separately using just one sample for each. 

First step is to trim and then do QC. Before this, you might want to run fastqc on the raw fasta files. 

```{bash}
#!/bin/bash
#SBATCH --job-name=bbduk
#SBATCH --account=rrg-ben
#SBATCH --array=0-9
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --output=bbduk.%J.%a.out
#SBATCH --error=bbduk.%J.%a.err
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load bbmap/39.06

in=/home/samp/projects/rrg-ben/for_Sam/muel/muel_fq/raw_fq
out=/home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed

declare -a forward=(${in}/*R1_001.fastq)
R1=${forward[${SLURM_ARRAY_TASK_ID}]}

adaptors=/home/samp/projects/rrg-ben/for_Sam/muel/ben_scripts/TruSeq2_and_3-PE-2.fa

base=$(basename ${R1} R1_001.fastq)

bbduk.sh threads=8 \
        in1=${in}/${base}R1_001.fastq \
        in2=${in}/${base}R2_001.fastq \
        out1=${out}/${base}R1_bbduk.fastq.gz \
        out2=${out}/${base}R2_bbduk.fastq.gz \
        ref=${adaptors} \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=36 ftl=1 maq=20

module load StdEnv/2023 fastqc/0.12.1

fastqc ${out}/${base}R1_bbduk.fastq.gz ${out}/${base}R2_bbduk.fastq.gz


```
This is my adaptor file:
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
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>PrefixPE/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA
>NEBNext_3p
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>NEBNext_5p
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Clontech_SMART_CDS_primer_II_A
AAGCAGTGGTATCAACGCAGAGTAC
>Clontech_SMART_CDS_primer_II_A_polyT
CAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```
Now you can actually run the assembly. This step is super time and memory intensive. (Edit: for me, this only took about 6.5 hours.)
```
#!/bin/sh
#SBATCH --job-name=trinity2
#SBATCH --cpus-per-task=32
#SBATCH --time=168:00:00
#SBATCH --mem=250G
#SBATCH --output=trinity2.%J.out
#SBATCH --error=trinity2.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2020
module load gcc/9.3.0 openmpi/4.0.3
module load trinity/2.14.0 samtools/1.17 jellyfish/2.3.0
module load salmon/1.4.0
module load python/3
module load scipy-stack/2023a

Trinity --seqType fq \
        --max_memory 200G \
        --CPU ${SLURM_CPUS_PER_TASK} \
        --full_cleanup \
        --min_kmer_cov 2 \
        --bflyCalculateCPU \
        --jaccard_clip \
        --left /home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed/X_muelleri_tad31_S11_L001_R1_bbduk.fastq.gz \
        --right /home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed/X_muelleri_tad31_S11_L001_R2_bbduk.fastq.gz \
        --output /home/samp/projects/rrg-ben/for_Sam/muel/new_improved_trinity2

#the --jaccard_clip flag is expensive but can help with dense genomes, such as those with many paralogs
#also, the order that you load modules is important... load compilers and such first... and don't forget to load python and scipy-stack so that you can use numpy!!!

```
When my assembly finished running, I had ~110000 unique transcripts. You can also check the gene_trans_map file to find out how many genes Trinity thinks you have (Trinity assembles isoforms independently for each gene). I had ~80000 genes. Keep in mind Trinity's identification of unique genes and isoforms is probably not very accurate. 80000 feels like a high estimate for the number of genes in X. muelleri: 
```
cut -f 1 trinity_tad31.fasta.gene_trans_map | sort | uniq | wc -l
```
Later, I can use the .gene_trans_map file as a tx2gene file with tximport in order to collapse isoforms to transcripts.   

It's important to run some kind of QC on the assembly to see how well it was assembled. Trinity offers a script for the ExN50 metric, which is a modification of the N50 metric. From a review by Raghavan et al: "The N50 is a simple metric which describes the sequence length at which half the nucleotides in the genome assembly are in sequences equal in or longer than this length [76]. The goal would then be to maximize the N50 value as this would indicate complete assembly of all genomic elements This is inappropriate for transcriptome assemblies as the objective is recovery of many (relatively) short full-length sequences, and not the construction of a few very long contigs ... . The ExN50 metric is a modification to the traditional N50 making it suitable for assessing transcriptome assemblies. Here, the N50 value is calculated only for the top X% of the cumulative expression levels. The length reported as corresponding to ExN50 is a ‘gene’ length obtained as the expression-weighted sum of the corresponding isoform lengths."   
   
In order to compute ExN50, you first have to quantify transcript abundance. You can do this with a Trinity script. The Trinity [documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification) shows multiple ways to get these estimates, but I chose kallisto since it's typically quite fast. 
```
#!/bin/sh
#SBATCH --job-name=trinity_abundance_est
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=250gb
#SBATCH --output=trinity_abundance_est.%J.out
#SBATCH --error=trinity_abundance_est.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2020 gcc/9.3.0 samtools/1.17 kallisto/0.46.1

/home/samp/projects/rrg-ben/for_Sam/trinity_ben/trinityrnaseq-v2.15.2/util/align_and_estimate_abundance.pl \
        --transcripts /home/samp/projects/rrg-ben/for_Sam/muel/trinity_tad31/trinity_tad31.fasta \
        --seqType fq \
        --left /home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed/X_muelleri_tad31_S11_L001_R1_bbduk.fastq.gz \
        --right /home/samp/projects/rrg-ben/for_Sam/muel/bbduk_trimmed/X_muelleri_tad31_S11_L001_R1_bbduk.fastq.gz \
        --est_method kallisto \
        --trinity_mode \
        --prep_reference \
        --output_dir /home/samp/projects/rrg-ben/for_Sam/muel/trinity_tad31/abundance_estimates

```
Once you have abundance estimates you can compute ExN50
```



```






Next, I want to assemble a transcriptome for a male individual. Once I have this assembly, I can use reciprocal best blast hits to collapse the transcriptomes together in order to find sex-shared and sex-specific transcripts. In order to assemble the male transcriptome, just choose the male sample with the most reads, and run the assembly script from above.  
   
Before doing BLAST reciprocal best hits, I filtered out just the longest isoform for each transcript: 
```
#!/bin/sh
#SBATCH --job-name=get_longest_transcript
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=100G
#SBATCH --output=get_longest_transcript.%J.out
#SBATCH --error=get_longest_transcript.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

/home/samp/projects/rrg-ben/for_Sam/trinity_ben/trinityrnaseq-v2.15.2/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/samp/projects/rrg-ben/for_Sam/muel/trinity_tad31/trinity_tad31.fasta > /home/samp/projects/rrg-ben/for_Sam/muel/trinity_tad31/trinity_tad31_longest_isoform.fasta
```
Now the reciprocal best hits. First make a blastable database for each transcriptome: 
```
makeblastdb -in trinity_tad31_longest_isoform.fasta -dbtype nucl -out trinity_tad31_longest_is
oform_blastable

```
Now that I have a single transcriptome with female-specific transcripts, male-specific transcripts and sex-shared transcripts, I want to run kallisto to quantify reads. Since I have a much smaller transcriptome, I should have much more statistical power: 
```
kallisto script will be here soon
```
