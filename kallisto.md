Since I'm regenerating kallisto counts using my newly trimmed reads and newly assembled transcriptome, I'm documenting how I use kallisto. 
The first step is to generate an index for my de novo transcriptome. The newest assembly of the transcriptome is located here: 
```
/home/samp/projects/rrg-ben/for_Sam/muel/sam_trinity_assembly/Trinity.fasta
```
This is how I generate the index: 
```
#!/bin/sh
#SBATCH --job-name=kallisto
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --output=kallisto.%J.out
#SBATCH --error=kallisto.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# run by passing an argument like this
# sbatch wip_kallisto_index path_to_transcriptome_assembly/transcriptomeAssembly.fasta

module load StdEnv/2020 intel/2020.1.217 kallisto/0.46.1

kallisto index -i ${1}.idx ${1}
```
The next step is to generate the quant files. This is my version of [Ben's script](https://github.com/evansbenj/2021_tad_ko_RNAseq/blob/75fa685bf62f5a7f5ebea7bf39b04f0cd6326fdf/kallisto.md#L4), adapted to run on all the samples in parallel using an array:
```
#!/bin/sh
#SBATCH --job-name=kallisto
#SBATCH --cpus-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=32gb
#SBATCH --output=kallisto.%J.out
#SBATCH --error=kallisto.%J.err
#SBATCH --account=def-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# run by passing an argument like this
# sbatch wip_kallisto_count ref_transcriptome_indexfile.idx folderwithfqfilez

module load StdEnv/2020 intel/2020.1.217 kallisto/0.46.1

index=${1}
dir=${2}
declare -a in=(${dir}/*R1.fq.gz)
fwd=${in[${SLURM_ARRAY_TASK_ID}]}
rev=${fwd/R1.fq.gz/R2.fq.gz}

if [[ -e ${fwd} && -e ${rev} ]]; then
        kallisto quant -b 100 -i ${ref} -o ${fwd/R1.fq.gz/kallisto_boot_out} ${fwd} ${rev}
fi
```
The final step is to compile all the quant files into a matrix: 
```

```
Next, you can use an R package like DESeq2 to analyze differential expression.
