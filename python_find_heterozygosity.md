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

#make sure the the number of individuals 



