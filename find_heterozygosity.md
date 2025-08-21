Annotated (by me) version of Ben's script, which is described as: "if a parent is not available (cal) this script can get positions that are homozygous in all individuals of one sex and heterozygous in at least (5) individuals of the other sex (separately for males and females)"
```
#!/usr/bin/env perl
use strict;
use warnings;
#"strict" forces you to declare variables using "my"
#"warnings" enables Perl's warning system 

# This program reads in a vcf file with genotypic information from
# a family and identifies positions that
# are homozygous in sons (mat sites because daughters are heteroz) and positions
# that are not homozygous only in daughters/mom (pat sites because sons are heteroz). 
# This is for dataset 
# where there is no information from one or both parent (e.g. the dad for calcaratus)
# module load StdEnv/2023 perl/5.36.1
# execute like this:
# ./Gets_matonly_positions_from_vcf_file_nopat.pl vcf 21111111111111111111000000000000000 matout patout

# epitrop RADseq: 111111111110000000000
# calcaratus RADseq: 11111111111111111111000000000000000

# where 21111111111111111111000000000000000 is whether the individuals are the excluded(2), a daughter(1), or a son(0)
# and matout and patout are the output files (with chromsoome name in them) 

my $vcf = $ARGV[0];
my $sexes = $ARGV[1];
my $outputfile = $ARGV[2];
my $outputfile2 = $ARGV[3];
#$ARGV[0] ect. are the built-in arrays for commandline arguments (since it is an array, the first argument is denoted[0]). The first argument is the vcf file you want to analyze, the second is a string of characters that represents the sexes (strings of 1s and 0s in the same order as the samples that denote male or female). The third and fourth arguments are the names you want for two different outfiles.

my @columns;
my @pat;
my @pat1;
my $x;
#these are local variables that you're going to use later

my @sexes = split("",$sexes);
#this is going to split the sexes variable (which is a string of numbers) into an array of single characters, each one representing the sex of a sample

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";
#"unless (...)" basically means "if the conditions in the round brackets aren't met, do x". In this case, if the outfile can't be opened, then print "I can't write to $outputfile" and $!, which denotes the error message associated with the last command, which was the attempt to open the outfile. Then the script exits.
#if the script WAS able to open or newly create the outfile, it prints "Creating output file $outputfile". Additionally, if the outfile already exists, the redirection operator ">" ensures that all the content in the file is deleted, so it can start out blank.

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2  $!\n\n";
	exit;
}
print "Creating output file: $outputfile2\n";
#same thing for the second output file


if ($vcf =~ /.gz$/) { 
	#open DATAINPUT, '<:gzip', $vcf or die "Could not read from $vcf: $!";
	open(DATAINPUT, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
}
#check if the vcf file ends with ".gz", indicating that it is gzipped. If so, open it using "gunzip" and pipe the decompressed output directly to the script. The "||" is the logical OR operator in Perl. It means that if the vcf cannot be decompressed, read and piped to the script, the script should "die" (exit with an error message) and print "Could not read from $vcf" with the associated "$!", which you should recall is the error message associated with the last command (which was the attempt to decompress the file)

else {
	open DATAINPUT, $vcf or die "Could not read from $vcf: $!";
}
#if the file does not have a .gz ending, assume it is not compressed and simply open it or "die" (exit with an error message)

my $number_of_samples=0;
my $switch1=0;
my $switch2=0;
my $num_daughters=0;
my $num_het_daughters=0;
my $num_sons=0;
my $num_het_sons=0;
#declare a bunch of local variables and start by setting them all to zero.

while ( my $line = <DATAINPUT>) {
	@columns=split("	",$line);
		#print $line,"\n";
		if($columns[0] =~ m/^#/){ # this is a commented line
			if($columns[0] eq '#CHROM'){ # this is the first line
				$number_of_samples = scalar(@columns)-9;
				print "Number of samples ",$number_of_samples,"\n";
			}
#the while loop continues processing the vcf file until it runs out of lines
#first, each line is split according to tabs and stores the parts in the "columns" array, this way the different fields are accessible each in their own array.
#columns[0] is the first field. Check if each line in the first field starts with a "#", which would indicate that they are headers. If one of these headers is "#CHROM", this means the script has found the last header line. You can use this last header line to calculate how many samples you have in your data. Count the scalar number of columns and subtract 9 (for the 9 default vcf columns), leaving you with the number of columns that are samples. Finally, print the number of samples
		}
		else{ # this is the genotype data #this "else" part of the loop processes lines that DON'T start with "#", meaning they are not headers. This loops through EACH variant site, which is represented by a single line in the vcf file. For each line, the script wants to figure out the genotype at this variant site of all sons and all daughters.
			$switch1=0; #tracks the genotyping status of daughters. Gets reset for each new variant site
			$switch2=0; #tracks the genotyping status of sons. Gets reset for each new variant site
			$num_daughters=0;
			$num_het_daughters=0;
			$num_sons=0;
			$num_het_sons=0; #these are all variables that i will start by setting to zero
			for ($x = 1 ; $x <= $number_of_samples; $x++ ){ # cycle through each sample, starting with sample 1. The sample number is "x", the "xth sample". "$x++" means that x increases by 1 each iteration of the loop. The loop continues so long as x is less than or equal to the number of samples, ensuring that all samples are genotyped. 
				if($sexes[$x-1] != 2){ # only consider samples that are included ("2" in the sexes string means to exclude a sample)
					@pat = split(":",$columns[$x+8]); # this is the whole genotype and other info from the vcf file. "$columns[$x+8]" gets the column for the current sample (the sample number x + 8 in order to skip the default vcf columns). The genotype info is colon-separated and gets stored in an array called @pat
					@pat1 = split(/[\|\/]/,$pat[0]); #this gets just the two alleles from the genotype string. This line split by either an "\" (phased allele) or an "|" (unphased allele). The alleles are stored in the @pat1 array
					# check if the daughters are homozygous or missing. This is the daughter part of the loop.
					if(($pat[0] ne './.')&&($pat[0] ne '.|.')&& 
					($pat[0] ne '.')&&($sexes[$x-1] == 1)){ # this is a daughter with a genotype. So this part of the loop only executes if the current sample is "1" in the sexes array (meaning female) and if the genotype is NOT missing
					# this is a daughter
					$num_daughters+=1; 
						if(($pat1[0] eq $pat1[1])&&($switch1 != 9)){ #if $pat1[0] is equal to $pat1[1], then the genotypes are the same and the daughter is homozygous at this position. If $switch1 != (not equal) to 9, then a heterozygous daughter has not been found yet. You want to find sites where ALL daughters are homozygous
							$switch1 = 1; #reset the switch to 1, meaning all daughters so far have been homozygous
							#print "goodpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch1 != 9)){ #if the alleles are not equal and switch1 does not equal 9, then this is the FIRST heterozgous daughter
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$switch1 = 9; #set switch1 to 9, indicating that you HAVE found a het daughter
							$num_het_daughters+=1; #increase the number of het daughters by one
						}	
						elsif(($pat1[0] ne $pat1[1])&&($switch1 == 9)){ #if the alleles are not the same, then this daughter is het for this position, but is not the first het daughter that has been found
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$num_het_daughters+=1;
						}	
					}
					# check if the sons are heterozygous or missing. This is the sons part of the loop. It functions the same as the daughter part, but it checks for the value of items in the sexes array to be "0" meaning male
					elsif(($pat[0] ne './.')&&($pat[0] ne '.|.')&&
					($pat[0] ne '.')&&($sexes[$x-1] == 0)){
						# this is a son
						$num_sons+=1;
						if(($pat1[0] eq $pat1[1])&&($switch2 != 9)){
							# this son is homoz 
							$switch2 = 1;
							#print "goodpats ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch2 != 9)){
							# this son is the first heteroz
							#print "badpats ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$switch2 = 9;
							$num_het_sons+=1;
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch2 == 9)){
							# this son is heteroz, but not the first one
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$num_het_sons+=1;
						}	
					}	
				}
				#print "\n";
			}
			#if(($switch1 == 1)&&($switch2 == 9)&&($num_sons == $num_het_sons)){
				# all the daughters with genotypes are homoz
				# all the sons with genotypes are heterozygous
			#if(($switch1 == 1)&&($switch2 == 9)){	
				# all the daughters with genotypes are homoz
				# at least one son with a genotypes is heterozygous
			if(($switch1 == 1)&&($switch2 == 9)&&($num_het_sons >=5)){		
				# all the daughters with genotypes are homoz
				# at least five sons with a genotypes are heterozygous
				print OUTFILE2 $columns[0],"\t",$columns[1],"\n";
				#outfile 2 contains the positions where all daughters are homozygous and at leas 5 sons are heterozygous
			}
			#elsif(($switch2 == 1)&&($switch1 == 9)&&($num_daughters == $num_het_daughters)){
				# all the sons with genotypes are homoz
				# all the daughters with genotypes are heterozygous
			#elsif(($switch2 == 1)&&($switch1 == 9)){
				# all the sonz with genotypes are homoz
				# at least one daughter with a genotype is heterozygous
			elsif(($switch2 == 1)&&($switch1 == 9)&&($num_het_daughters >=5)){				
				# all the sonz with genotypes are homoz
				# at least five daughters with a genotype are heterozygous
				print OUTFILE $columns[0],"\t",$columns[1],"\n";
				#outfile 1 contains the positions where all the sons are homozygous and at leas five daughters are heterozygous. This is closest to the conditions I want for muelleri, although I probably want ALL daughters with genotypes to be heterozygous and all sons with genotypes to be homozygous.
			}
			
		} # end else
} # end while

```
close DATAINPUT;
close OUTFILE;
close OUTFILE2;
