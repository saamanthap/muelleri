I'm using RepeatMasker to identify conserved repeats in the 5kb region upstream of ddx4 (this is the vasa homologue that I am using to make a reporter construct). In order to annotate species-specific repeats, you can download custom databases from dfam.org: https://www.dfam.org/releases/current/families/FamDB/ You must install the root first, then you can add any of the other partitians, which are separated by taxa. I'm downloading partition 12 which contains repeats for frogs and toads. RepeatMasker is already installed on Nibi, but not all partitions are present.

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
Then when you run the RepeatMasker command itself, you need to specify what database you want to use. (My input file here was a fasta file with entries for the 5'UTR and 5kb upstream region for ddx4 in tropicalis and laevis L and S subgenomes.)

```
blah blah blah RepeatMasker script
```

