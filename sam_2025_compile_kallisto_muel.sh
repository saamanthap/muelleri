Modified version of 2025_compile_kallisto_muel.sh that automatically locates abundance.tsv files within a given directory.

#!/bin/bash
#SBATCH --job-name=compile_kallisto_counts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=kallisto.%J.out
#SBATCH --error=kallisto.%J.err
#SBATCH --account=rrg-ben

module load StdEnv/2020 trinity/2.14.0


# run like this
# copy_sbatch 2025_compile_kallisto_muel.sh directory_name

# if the number of arguments passed to the script is not exactly one, the program prints an error message and exits
if [ "$#" -ne 1 ]
then
        echo "Usage: $0 <directory>"
        exit 1
fi

# the value of the first command line argument is substituted into the $1 string and this value is stored in the target_directory variable
target_directory="$1"

# this will check whether the directory passed to the string is a valid directory
if [ ! -d "$target_directory" ]
then
        echo "Error: $target_directory is not a valid directory."
        exit 1
fi

# this defines the name of the file i want to read from each subdirectory
filename="abundance.tsv"

#define a string that I'm going to concatenate
commando="./abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix muel_kallisto_countz_ --gene_trans_map none --name_sample_by_basedir"

# find all the subdirectories inside the given directory (only look one level deep)(this is accomplished using process redirection at the end of the loop)
# for every subdirectory that the program is able to read, construct the path to each abundance.tsv file and define the path as the filepath variable
while read -r subdirectory
do
        filepath="$subdirectory/$filename"

        if [ -f "$filepath" ]
        then
                echo "Reading file: $filepath"
                commando="$commando $filepath"
        else
                echo "Error: $filename not found in $subdirectory"
        fi
done < <(find "$target_directory" -type d -mindepth 1 -maxdepth 1)
# the while loop ends here

# execute the string that has now been built
eval "$commando"

# print a confirmation that the entire script was run successfully
echo "Yippee!"
exit 0
