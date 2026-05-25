I'm using RepeatMasker to identify conserved repeats in the 5kb region upstream of ddx4 (this is the vasa homologue that I am using to make a reporter construct). In order to annotate species-specific repeats, you can download custom databases from dfam.org: https://www.dfam.org/releases/current/families/FamDB/ You must install the root first, then you can add any of the other partitians, which are separated by taxa. I'm downloading partition 12 which contains repeats for frogs and toads (as well as a bunch of other taxa including weird jawless fish). RepeatMasker is already installed on Nibi, but not all partitions are present.

I ran the download as a slurm submission script since some paritions are very large: 

```
#!/bin/bash
#SBATCH --job-name=download_dbfam
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=24gb
#SBATCH --output=download_dbfam.%J.out
#SBATCH --error=download_dbfam.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL


wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz
gunzip dfam39_full.0.h5.gz

wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.12.h5.gz
gunzip dfam39_full.12.h5.gz

```
In order to make use of these databases, you need to extract your desired repeats as a fasta file. I'm replicating what was done in this paper, so I extracted anura (all frogs and toads) repeats only. I found this step to be time-consuming, so I submitted it to the queue as a job. The only contents of the directory specified after -i should be the root and other partitions you want to use.

```
#!/bin/bash
#SBATCH --job-name=extract_fasta_famdb
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=24G
#SBATCH --output=extract_fasta_famdb.%J.out
#SBATCH --error=extract_fasta_famdb.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2023 gcc/12.3 repeatmasker/4.2.1

famdb.py -i /home/samp/projects/rrg-ben/for_Sam/reporter_project/famdb families -f fasta_name -d anura > ${1}
```

Then when you run the RepeatMasker command itself, you need to specify what database you want to use. (My input file here was a fasta file with entries for the 5'UTR and 5kb upstream region for ddx4 in tropicalis and laevis L and S subgenomes.)

```
blah blah blah RepeatMasker script
```

