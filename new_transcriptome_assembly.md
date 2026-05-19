# New and improved transcriptome assembly!
I am trying to assemble a new transcriptome where I have less independently assembled alleles. This will increase my statistical power when it comes to the differential expression analysis. The idea here is that I will trim more thoroughly (using bbduk instead of trimmomatic and cutadapt), and then assemble a male and female transcriptome separately using just one sample for each. 

First step is to trim and then do QC. Before this, you might want to run fastqc on the raw fasta files. 

```{bash}
#!/bin/bash
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

base=`basename ${R1}R1_001.fastq`

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
```
