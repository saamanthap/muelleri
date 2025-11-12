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
In order to use these scripts, I have to convert my .vcf files to .geno files (basically a tab-delimited file that contains minimal information). You can use the script parseVCFs.py in the VCF_processing directory to do this. (It can convert multiple vcfs at once.) Run it like this:
```
python parseVCFs.py -i input.vcf.gz -1 input2.vcf.gz | bgzip > output.geno.gz
```
I'm using the popgenWindows.py script to estimate Fst for fixed windows on the chromosome. You can run this:
```
