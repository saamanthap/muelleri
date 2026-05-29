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
Once you have abundance estimates you can compute ExN50. You can do this at a gene or transcript level (according to Trinity documentation, gene is more useful). Below I show the gene-level analysis, but you can compute ExN50 for individual transcripts by changing the word "gene" to transcript.
```
#!/bin/sh
#SBATCH --job-name=trinity_ExN50
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=24gb
#SBATCH --output=trinity_ExN50.%J.out
#SBATCH --error=trinity_ExN50.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL


/home/samp/projects/rrg-ben/for_Sam/trinity_ben/trinityrnaseq-v2.15.2/util/misc/contig_ExN50_statistic.pl \
        /home/samp/projects/rrg-ben/for_Sam/muel/trinity_tad31/abundance_estimates/abundance.tsv \
        /home/samp/projects/rrg-ben/for_Sam/muel/trinity_tad31/trinity_tad31.fasta \
        gene | tee tad31_ExN50.gene.stats

#you can do this at the transcript level or gene level... just use the command "gene" or "transcript" before the pipe to specify the target

```
Here is an example of the gene-level output. A good way to interpret this is by plotting Ex against the ExN50 value. For a well-assembled transcriptome, N50 should peak around Ex90. Since mine crashes early, only the most highly expressed transcripts are assembled well, and then rest are highly fragmented:

```
Ex      ExN50   num_genes
1       2471    12
2       3947    30
3       4626    55
4       4671    85
5       4640    120
6       4740    158
7       4600    201
8       4462    246
9       4450    295
10      4218    345
11      4160    398
12      4045    453
13      3971    511
14      3916    572
15      3883    635
16      3826    701
17      3757    769
18      3731    840
19      3654    914
20      3617    989
21      3595    1067
22      3557    1149
23      3522    1232
24      3501    1318
25      3462    1407
26      3432    1499
27      3390    1593
28      3360    1691
29      3338    1791
30      3320    1895
31      3291    2002
32      3244    2112
33      3198    2226
34      3174    2343
35      3133    2464
36      3106    2589
37      3087    2719
38      3064    2853
39      3047    2990
40      3026    3133
41      3006    3280
42      2980    3434
43      2947    3593
44      2924    3758
45      2892    3929
46      2863    4107
47      2835    4291
48      2801    4483
49      2779    4683
50      2751    4891
51      2731    5108
52      2710    5335
53      2687    5572
54      2672    5821
55      2652    6081
56      2636    6354
57      2622    6641
58      2610    6942
59      2596    7260
60      2584    7595
61      2568    7949
62      2546    8324
63      2509    8720
64      2468    9138
65      2427    9584
66      2390    10059
67      2350    10566
68      2308    11107
69      2268    11687
70      2226    12312
71      2187    12987
72      2145    13712
73      2102    14493
74      2062    15337
75      2024    16253
76      1985    17247
77      1943    18329
78      1901    19503
79      1861    20786
80      1821    22185
81      1782    23715
82      1736    25385
83      1694    27207
84      1659    29190
85      1615    31344
86      1568    33679
87      1523    36203
88      1479    38925
89      1434    41846
90      1393    44965
91      1353    48285
92      1310    51814
93      1269    55527
94      1230    59410
95      1188    63468
96      1148    67711
97      1107    72153
98      1067    76815
100     1026.19016255607 
```
I want to figure out if this new assembly is better or worse than my previous assembly, which was assembled from all ten samples. It's possible that having more samples creates more "depth" which could help me assemble transcripts more completely. I already have an abundance file for this transcriptome (I had to estimate counts for all samples using kallisto and then combine the counts into a single file with a column for each sample).

```
#!/bin/sh
#SBATCH --job-name=trinity_ExN50
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=24gb
#SBATCH --output=trinity_ExN50.%J.out
#SBATCH --error=trinity_ExN50.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL


/home/samp/projects/rrg-ben/for_Sam/trinity_ben/trinityrnaseq-v2.15.2/util/misc/contig_ExN50_statistic.pl \
        /home/samp/projects/rrg-ben/for_Sam/muel/kallisto_fully_trimmed/muel_kallisto_countz_fully_trimmed.isoform.TMM.EXPR.matrix \
        /home/samp/projects/rrg-ben/for_Sam/muel/sam_trinity_assembly/Trinity.fasta \
        gene | tee trinity_ExN50.gene.stats

#you can do this at the transcript level or gene level... just use the command "gene" or "transcript" before the pipe to specify the target

```
It turns out this transcriptome is *much* less fragmented, despite having more than 560000 transcripts. The ExN50 value peaks at around Ex93. Apparently these many transcripts are not fragmented, they are just independantly assembled alleles or isoforms. 



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
