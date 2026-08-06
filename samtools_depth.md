Compute depth for all positions using samtools depth:
```
#!/bin/bash
#SBATCH --job-name=samtools_depth
#SBATCH --cpus-per-task=1
#SBATCH --array=0-31
#SBATCH --time=100:00:00
#SBATCH --mem=48gb
#SBATCH --output=samtools_depth.%J.%a.out
#SBATCH --error=samtools_depth.%J.%a.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2023 samtools/1.22.1

bam_dir=${1}
ending=${2}

declare -a bams=(${bam_dir}/*"${ending}")

current_bam=${bams[${SLURM_ARRAY_TASK_ID}]}
base=$(basename ${current_bam} ${ending})

samtools depth -aa -o ${base}_samtoolsdepth.out ${current_bam}

```
Merge all depth files in a memory-efficient manner. This only works if you used flag -aa to compute depth for ABSOLUTELY ALL positions. The file this produces tends to be HUGE (>250gb). You should zip all output files to save space.
```
#!/bin/bash
#SBATCH --job-name=merge_depth
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=100gb
#SBATCH --output=merge_depth.%j.out
#SBATCH --error=merge_depth.%j.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

# track problems in case of failure
set -euo pipefail

# outfile name
out=${1}
files=( *samtoolsdepth.out )

# set up temp dirs
tmpdir=$(mktemp -d -p ${SLURM_TMPDIR:-/tmp})
trap 'rm -rf "$tmpdir"' EXIT

# build the header
header="CHROM\tPOS"
for f in ${files[@]}; do
        sname=$(basename $f | sed 's/\.samtoolsdepth\.out$//')
        header="${header}\t${sname}"
done

# write the header to the output file
echo -e "$header" > ${out}

# extract depth column from each file to temp storage
cut_files=()
for ((i=1; i<${#files[@]}; i++)); do
        tmp_cut="${tmpdir}/col_${i}.txt"
        cut -f 3 "${files[$i]}" > "$tmp_cut"
        cut_files+=( "$tmp_cut" )
done

# paste first file and add the depth extracted from remaining files
paste "${files[0]}" "${cut_files[@]}" >> "${out}"
```
Now I need to extract positions with zero 
