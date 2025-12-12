## Installing pblat-cluster
Download the source code to your desired directory:
```
git clone https://github.com/icebert/pblat-cluster
```
Enter the source code directory and use the command "make"
```
cd pblat cluster
make
```
If the make command doesn't work, you may have to modify the Makefile. Take this line: 
```
$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) jkweb.a $(MPILIBS)
```
Change it to this:
```
$(CC) -Wl,--allow-multiple-definition $(LDFLAGS) -o $@ $(OBJS) $(LIBS) jkweb.a $(MPILIBS)
```
Then clean and recompile:
```
clean make
make
```
## Using pblat-cluster
pblat-cluster is a version of blat for mpi-compatible clusters (like Nibi!). Basically, it adds parallelization to blat. I used the following script to align my de novo transcriptome assembly to the borealis genome:
```
#!/bin/bash
#SBATCH --job-name=pblat-cluster
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --time=23:00:00
#SBATCH --mem=200gb
#SBATCH --output=pblat-cluster.%J.out
#SBATCH --error=pblat-cluster.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN,END,FAIL

module load StdEnv/2023 openmpi/4.1.5

database=${1}
query=${2}
output=${3}
pblat=/home/samp/pblat-cluster/pblat-cluster


mpirun ${pblat} -q=rnax -t=dnax -out=blast8 ${database} ${query} ${output}
```
Note: -q=rnax -t=dnax are the suggested flags for aligning full-length mRNA against a genome from another species.
