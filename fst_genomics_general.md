I'm using the [genomics_general](https://github.com/simonhmartin/genomics_general) collection of scripts to estimate Fst values for sites in my samples. First clone the repository to my home directory on nibi: 
```
cd #go to home directory
git clone https://github.com/simonhmartin/genomics_general.git
```
Install numpy in your python environment if you don't have it already
```
module load python
source ~/my_env/bin/activate
pip install numpy
```
I also had to make a change to the genomics.py source code because the np.NaN attribute has become np.nan in recent version of numpy. I did this very quickly in vi: 
```
vi genomics.py
:%s/np.NaN/np.nan/g
```
In order to use these scripts, I have to convert my .vcf files to .geno files (basically a tab-delimited file that contains minimal information). You can use the script parseVCF.py in the VCF_processing directory to do this. You can run it on every vcf at once directly from the commandline (this could take a few minutes):
```
for i in $(ls *_out.vcf) ; do python ~/genomics_general/VCF_processing/parseVCF.py -i $i > ${i}.geno ; done
```
I'm using the popgenWindows.py script to estimate Fst for fixed windows on the chromosome. This is the command I'm toying with right now:
```
python ~/genomics_general/popgenWindows.py -g DBI_combineGVCF_Chr1L_out.vcf.geno -o DBI_combineGVCF_Chr1L_out.vcf.geno.out --windType coordinate -f phased -w 1000000 -m 100 -p F -p M --popsFile pops_file.txt
```
Note: the script seems to struggle to compute fst unless the -minSites flag (-m) is set very low... this is probably because there are quite a few windows with low numbers of sites...  

Here is a simple loop that allows me to process all the .geno files at once: 
```
for i in $(ls *_out.vcf.geno) ; do python ~/genomics_general/popgenWindows.py -g DBI_combineGVCF_Chr1L_out.vcf.geno -o ${i}.out --windType coordinate -f phased -w 1000000 -m 100 -p F -p M --popsFile pops_file.text ; done 

```





