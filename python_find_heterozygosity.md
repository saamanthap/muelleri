Here is my version of Ben's script to find female-specific heterozygosity from a vcf file:
```
#!/usr/bin/env python3  #this is where the correct shebang will go
import sys #this library lets you use commandline arguments
import pysam #this library is for handling vcf files (so that i don't have to do the manual splitting that Ben used in his script)

def main():
  vcf_file = sys.argv[1] #this is the first commandline argument, giving the vcf file
  sex_string = sys.argv[2] #this is the second commandline argument, giving the string that denotes sex: 0010000110
  output_file = sys.argv[3] #this is the third commandline argument, giving the name of the first output file
  output_file2 = sys.argv[4] #this is the fourth commandline argument, giving the name of the second output file

  #this will take the sex_string and create a dictionary called "sexes_map".
  #i take the sex_string variable and use "enumerate" to give the index and the value of each character in the string.
  #int(s) is what builds the dictionary. "i" becomes the key and int(s) is the value (meaning the character from the string)
  sexes_map = {i: int(s) for i, s in enumerate(sex_string)}

  #use "with" to open files and automatically close them when done
  #give the files "handles" (basically nicknames), out_f and out_f2
  #pysam handles gzipped files automatically
  with open(output_file, 'w') as out_f, open(output_file2, 'w') as out_f2:
    #"try" is the beginning of a "try...except" block that is watched for errors. If the code in the "try" fails, then Python throws a "FileNotFoundError" exception that is caught by the "except" part of the code, which is specifically meant to listen for and handle the error.
    try:
      vcf_in = pysam.VariantFile(vcf_file, 'r')
    except FileNotFoundError:
      sys.stderr.write(f"Error: Could not open VCF file {vcf_file}\n")
      sys.exit(1)
      



  pass

if __name__ == "__main__":
  main()



```
