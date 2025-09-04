Here is my python script: 

```
#first i'm going to import useful libraries
import sys #allows me to exit on errors
import argparse #allows me to use commandline arguments
import pysam #allows me to interpret vcf files directly (rather than relying on converting to a tab-delimited format first)
from collections import Counter #import the "Counter" class from the "collections" module
#collections is a module that provides specialized datatypes that are ideal for specific tasks
#Counter is a subclass of dict that stores elements as dictionary keys and their counts as dictionary values - using a Counter is an efficient way to make a "tally" of elements

def main():
#i start by setting up the parser object which will control all the commandline arguments. i add a help message that can be printed out if the script is run with the --help flag
  parser = argparse.ArgumentParser(
    description='Find sex-linked variants from a VCF file.',
    formatter_class=argparse.RawTextHelpFormatter
)

#now i'm going to actually define all the arguments i want to use. the first one is the vcf file. i included a help message reminding me what file formats are acceptable
parser.add_argument('vcf_file', help='Input VCF file (.vcf or .vcf.gz).')
parser.add_argument('sex_string', help='String identifying sex of each individual: 0=male, 1=female, 2=ignore.\n'
  'The order must match the order of samples in the VCF file.'
)
parser.add_argument('output_file', help='Output file for sites of interest.')
#the float data type means a number that contains a decimal. in this case the proportion is a decimal between 0-1
parser.add_argument(
  'proportion',
  type=float,
  help='Proportion of individuals within a sex that must gave a genotype AND show the specified pattern.'
#the next line actually runs the parsing process. "parser.parse_args()" will return an object that contains all the parsed values, which will be stored in the "args" variable
args = parser.parse_args()

#in order to access the parsed arguments, i'm going to store the values that are returned in variables with descriptive names
vcf_file = args.vcf_file
sex_string = args.sex_string
output_file = args.output_file
proportion = args.proportion

#now i need to parse and validate the input (basically read the input and make sure the input is appropriate)
#first check if the proportion is between 0-1. i have to write 0 as 0.0 because proportion is a float.
#sys.stderr sends the error message to standard error, in order to keep it separate from the script's regular output
if not 0.0 <= proportion <= 1.0:
  print("Error: Proportion argument must be between 0 and 1.", file=sys.stderr)
  sys.exit(1)

#try to convert each of the characters in the sex string into integers, which are stored in a list. if that doesn't work, print an error message and exit. "ValueError" watches for values that are not integers
try:
  sexes = [int(s) for s in sex_string]
except ValueError:
  print("Error: Sex string must only contain 0, 1, or 2.", file=sys.stderr)
  sys.exit(1)

#try to open the vcf file using pysam. use pysam.VariantFile to open the vcf
try:
  vcf_in = pysam.VariantFile(vcf_file)
#IOError watches for errors that happen when the file doesn't exist, you don't have permission to read it, and other file-related errors. "as e" saves the specific error message so that you can print it out
except IOError as e:
  print(f"Error: Could not open VCF file {vcf_file}. {e}", file=sys.stderr)
  sys.exit(1)

#make sure the the number of samples matches the length of the sex string. recall that "sexes" is the sex string formatted as a list of integers. vcf_in.header.samples effectively represents the samples within the header of the "vcf_in" file
if len(sexes) != len(list(vcf_in.header.samples)):
  print("Error: Number of individuals in sex string does not match VCF samples.", file=sys.stderr)
  sys.exit(1)

#build indices that contain the male and female samples. basically, i'm building new lists of pairs (index, value) 
male_indices = [i for i, s in enumerate(sexes) if s == 0]
female_indices = [i for i, s in enumerate(sexes) if s ==1]

#count the males and females
num_males = len(male_indices)
num_females = len(female_indices)
print(f"This includes {num_females} female(s) and {num_males} male(s).")

#set up the output file. the 'w' in open() is write mode, meaning that the file is created if it doesn't exist, and if it does exist, it is overwritten
try:
  outfile = open(output_file, 'w')
  outfile.write("CHR\tPS\tTYPE\tCATEGORY\tn_FEM\tn_MALE\n")
except IOError as e:
  print(f"Error: Could not open output file {output_file} for writing. {e}", file=sys.stderr)
  sys.exit(1)

#now is the time to apply the main logic of the script! i will iterate through the vcf file
#"record" is a type of object specific to vcf files that has special properties for the chromosome, position and genotype information for every sample
#i will also initialize empty lists male_alleles and female_alleles. these lists get reset with each position
for record in vcf_in:
  male_alleles = []
  female_alleles = []

for sample_index in male_indices:
  sample_name = vcf_in.header.samples[sample_index] #use the index to get the name of a particular sample from the header
  alleles = record.samples[sample_name].alleles #use the sample name to get the allele for a particular sample at a particular position
  if None not in alleles:
    male_alleles.extend(alleles) #in the pysam library, missing genotypes are represented by "None". if the genotype is not missing, add the allele to the male_alleles list. The list will collect all the alleles for all the males for a given position, and will reset for the next position.

#do the same for female alleles
for sample_index in female_indices:
  sample_name = vcf_in.header.samples[sample_index]
  alleles = record.samples[sample_name].alleles
  if None not in alleles:
    female_alleles.extend(alleles)

#each genotyped individual has two alleles, so to get the number of genotyped individuals, divide the number of alleles by 2
#basically, this block is going to check if the number of genotyped males and females is 0 (meaning all individuals have missing genotypes). if true, the rest of the code is skipped for the position, since there is no point executing on empty lists
num_males_genotyped = len(male_alleles) / 2
num_females_genotyped = len(female_alleles) / 2
if num_males_genotyped == 0 and num_females_genotyped == 0:
  continue

#this is the part where i use counter! this replaces the way ben used "uniq" in the original perl script
#i'm going to create Counter objects that count the number of unique alleles that appear in the male_alleles and female_alleles lists
male_allele_counts = Counter(male_alleles)
female_allele_counts = Counter(female_alleles)

is_male_hom = len(male_allele_counts) == 1 #only 1 unique allele was found, so all males are homozygous at this position
is_male_het = len(male_allele_counts) > 1 #more than 1 unique allele was found, so at least one male is heterozygous
is_female_hom = len(female_allele_counts) == 1
is_female_het = len(female_allele_counts) > 1

#Now i want to check various scenarios that i am interested in, and save the positions where these scenarios occur
#first scenario: males are homozygous and females are heterozygous (what i would likely expect for a muelleri sex-determining region) called Category 1
if is_male_hom and is_female_het and num_males_genotyped > 0: #making sure that males are homozygous, females are heterozygous, and at least 1 male was genotyped. if so, check proportion of divergence

#here, i can get the single allele carried by the homozygous males. i will simply take the first item stored in the Counter object male_allele_counts's keys
male_base = list(male_allele_counts.keys())[0]
#count how many female alleles are different from the single male allele. use [male_base] to subset the female_allele_counts, essentially counting how many times the male allele appears and subtracting from the total number of female alleles
diverged_alleles_count = len(female_alleles) - female_allele_counts[male_base]

#i specify the proportion of divergence in a commandline argument when running the script. here, i check if the number of diverged alleles is greater than the proportion times the number of female alleles. if proportion were 0.5, this would check that at least half of the female alleles diverge from the male allele (meaning that all females are heterozygous!) <- not sure if this is a good explanation?
if diverged_alleles_count > proportion * len(female_alleles):
  outfile.write(
    f"{record.chrom}\t{record.pos}\tSex_specific_SNP\t1\t"
    f"{num_females_genotyped}\t{num_males_genotyped}\n"
)

#check for the extreme case where all females are heterozygous
#loop through female_alleles in increments of two and add 1 for each pair of alleles IF they are different (meaning the individual is heterozygous)
num_het_females = sum(1 for i in range(0, len(female_alleles), 2)
                        if female_alleles[i] != female_alleles[i+1])
if num_het_females == num_females_genotyped: #checks if ALL genotyped females are het and then writes the position to the outfile
  outfile.write(
    f"{record.chrom}\t{record.pos}\tSex_specific_heterozygosity\t1\t"
    f"{num_females_genotyped}\t{num_males_genotyped}\n"
  )

#Category -1 - female homozgosity and male heterozygosity
elif is_female_hom and is_male_het and num_females_genotyped > 0:
  female_base = list(female_allele_counts.keys())[0]
  diverged_alleles_count = len(male_alleles) - male_allele_counts[female_base]

  if diverged_alleles_count > proportion * len(male_alleles):
    outfile.write(
      f"{record.chrom}\t{record.pos}\tSex_specific_SNP\t-1\t"
      f"{num_females_genotyped}\t{num_males_genotyped}\n"
    )

#extreme case - all males are heterozygous
num_het_males = sum(1 for i in range(0, len(male_alleles), 2)
                      if male_alleles[i] != male_alleles[i+1])
if num_het_males == num_males_genotyped:
  outfile.write(
    f"{record.chrom}\t{record.pos}\tSex_specific_SNP\t-1\t"
    f"{num_females_genotyped}\t{num_males_genotyped}\n"
  )

#category 1 - fixed divergence (meaning both males and females are homozygous, but for different alleles)
elif is_male_hom and is_female_hom and male_allele_counts != female_allele_counts:
  outfile.write(
    f"{record.chrom}\t{record.pos}\tFixed_divergence\t1\t"
    f"{num_females_genotyped}\t{num_males_genotyped}\n"
  )

#category 1 - female specific nucleotide (meaning males have no data!)
elif num_males_genotyped == 0 and num_females_genotyped > 0:
  if num_females_genotyped > proportion * num_females:
    outfile.write(
      f"{record.chrom}\t{record.pos}\tSex_specific_nucleotide\t1\t"
      f"{num_females_genotyped}\t{num_males_genotyped}\n"
    )

#category -1 - male specific nucleotide (meaning females have no data)
elif num_females_genotyped == 0 and num_males_genotyped > 0:
  if num_males_genotyped > proportion * num_males:
    outfile.write(
      f"{record.chrom}\t{record.pos}\tSex_specific_nucleotide\t-1\t"
      f"{num_females_genotyped}\t{num_males_genotyped}\n"
    )

#now close the files
outfile.close()
vcf_in.close()

if __name__ == "__main__":
  main()






