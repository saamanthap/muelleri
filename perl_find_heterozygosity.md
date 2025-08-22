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
			
		} # end else (this is the else loop that continued as long as the script could read NON-header line)
} # end while (this is the while loop that continued so long as there were lines of the file to read)

close DATAINPUT;
close OUTFILE;
close OUTFILE2;

```
Another of Ben's scripts, which "screens for sex-specific hets". Here is what the tab-delimited file looks like (for my data, even though the script is designed for Ben's data): 
```
#CHROM  POS     REF     X_muelleri_tad31_S11_L001__trim_sorted.bam      X_muelleri_tad32_S12_L001__trim_sorted.bam  X_muelleri_tad33_S13_L001__trim_sorted.bam      X_muelleri_tad34_S14_L001__trim_sorted.bam  X_muelleri_tad35_S15_L001__trim_sorted.bam      X_muelleri_tad36_S16_L001__trim_sorted.bam          X_muelleri_tad37_S17_L001__trim_sorted.bam      X_muelleri_tad38_S18_L001__trim_sorted.bam          X_muelleri_tad39_S19_L001__trim_sorted.bam      X_muelleri_tad42_S20_L001__trim_sorted.bam
Chr4L   24979   T       ./.     ./.     ./.     G/G     ./.     ./.     ./.     ./.     ./.     G/G
Chr4L   40872   G       ./.     A/A     ./.     ./.     A/A     ./.     ./.     ./.     ./.     ./.
Chr4L   40877   T       ./.     G/G     ./.     ./.     G/G     ./.     ./.     ./.     ./.     ./.

```
The script (annotated by me)
```

#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
#  Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for sites that are homozygous
#  in one sex for one SNP and at least partially heterozygous in the other sex

# to execute type Parse_tab.pl inputfile.tab 1111100110000111100011100110010100002200 interesting_sites.out proportion
# where 1111100110000111100011100110010100002200 refers to whether or not each individual in the ingroup 
# in the vcf file is (0) male, (1) female, and or (2) skipped

# proportion is the proportion of genotyped alleles in the heterogametic sex that are required to be
# different from the homogametic sex in order for the position to be reported.  This is a way to reduce reporting
# of low frequency polymorphisms (which are unlikely to be sex-linked but likely to have one sex all homozygous).
# the proportion parameter should be less than or equal to 0.5 

# if it is 0.5, this means all females are heterozygous and all males are homozygous (for positions with only 2 variants)

# we will also use this proportion to be a requirement for male-specific or female-specific SNPs, meaning at least
# this proportion of the individuals within each sex is required to have a genotype.

# het_sites.out sex_specific_sites.out diverged_sites.out are the output files that have the positions and chr of interesting sites

# example for clivii
# perl Parse_tab.pl clivii_unfiltered_removed_allchrs.vcf.tab 111111111111111111111111110000000000000000000 interesting_sites.out 0.35
# include only Eritrea:
# 222222221111111111111112222222222200000000222

# exclude Eritrea:
# 111111112222222222222221110000000022222222000

# Example for XB_WGS
# perl Parse_tab.pl XB_WGS_not_filtered_allchrs.vcf.gz.tab 100110011101010000102222 interesting_sites.out 0.5

my $inputfile = $ARGV[0]; #this is the tab-delimited text file
my $input2 = $ARGV[1]; #this is the string that identifies sex
my $outputfile1 = $ARGV[2]; #this is the name of the outfile
my $proportion = $ARGV[3]; #this is the proportion (a decimal) of individuals within each sex that must have a genotype

print "hello ",$proportion,"\n";

unless (open DATAINPUT, $inputfile) { #print an error message UNLESS the file is sucessfully opened
	print "Can not find the input file.\n";
	exit;
} 

unless (open(OUTFILE1, ">$outputfile1"))  { #print an error message UNLESS you are able to open and redirect data to the outfile
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n"; #if nothing goes wrong, print the name of the output file 
print OUTFILE1 "CHR\tPOS\tTYPE\tCATEGORY\tn_FEMs\tn_MALS\n"; #print the outfile with the given info in a tab-delimited format




my @sexes = split("",$ARGV[1]); #split the sex string into characters, each representing the sex of an individual

my @males=(); #open empty string
my @females=(); #open empty string
my @temp;
my @unique_male_nucleotides;
my @unique_female_nucleotides;
my $y; #represents an "index number" (sort of) of which individual i'm currently working on
my $x; #represents the number of female individuals
my $counter=0; #declare variable, set to zero initially
my $diverged=0;
my $number_of_male_individuals_genotyped=0;
my $number_of_female_individuals_genotyped=0;

for ($y = 0 ; $y <= $#sexes ; $y++ ) { #do this loop for the zeroth individual, while the index number of the individual is less than or equal to the last index of the sex array. increase y by one each iteration.
	if($sexes[$y] == 0){ #if the sex of the current individual is male, add 1 to the number of male individuals that have been genotyped so far
		$number_of_male_individuals_genotyped +=1; 
	}	
}	
for ($y = 0 ; $y <= $#sexes ; $y++ ) { #do this loop for the zeroth individual, while the index number of the individual is less than or equal to the last index of the sex array
	if($sexes[$y] == 1){ #if the individual is female, add 1 to the number of female individuals that have been genotyped so far
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s) and  ", $number_of_male_individuals_genotyped," males\n"; #print the number of females and the number of males

while ( my $line = <DATAINPUT>) { #this loop is what actually processes the file line by line. The loop continues WHILE a line is successfully read from the file. The loop terminates when it runs out of new lines
	chomp($line); #chomp removes trailing newline characters from a string to make sure it only contains the info i want
	@temp=split /[\t\/]/,$line; #split the line according to tab-delimited fields and backslashes. Store the info in the temp string. This creates ((# of individuals)x2 +3) number of fields. 
	if($temp[0] ne '#CHROM'){ # "#CHROM" is part of the header, so this condition allows you to execute code only for the non-header lines. Specifically, this checks if the first column is not equal to #CHROM
		if($#temp ne (($#sexes+1)*2)+2){ #this checks the the largest index number of temp to make sure the amount of sample columns matches the number of characters in the sex string. if not, print an error message
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n"; 
		}

		# parse the bases in all genotypes in each sex
		@males=();
		@females=();
		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) { #y=3 is true (which represents the ref genotype column), the current column number is less than or equal to the total column number, y=y+2 means move over two columns (second genotype in the first individual's genotype column) 
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")){ #if the first allele is not missing and the second allele is not missing
				if($sexes[$counter] == 0){ #if counter variable is set to zero (male individual)
						push(@males, $temp[$y]); #push adds elements to the end of an array. in this case, it adds the value of the first allele for the current individual
						push(@males, $temp[$y+1]); #adds the value of the second allele for the current individual
				}
				elsif($sexes[$counter] == 1){ #if the current individual is female
					push(@females, $temp[$y]); #add the first allele to the females array
					push(@females, $temp[$y+1]); #add the second allele
				}	
			}
			$counter+=1; #increase the value of the counter variable by 1 (to indicate that you have finished loading all the non-missing genotypes into the @males and @females arrays)
		}	
		# OK I should have all the bases loaded for non-missing genotypes for each male and each female
		
		@unique_male_nucleotides = uniq @males; #extract only the unique elements of the males array, meaning only look at the heterozygous positions. Recall that this array contains alleles only (just letters). After applying "uniq", heterozygous positions should have at least 2 items in the array (because there were two different alleles)
		@unique_female_nucleotides = uniq @females; #extract only the unique elements of the females array
		#print @females," ",@males,"\n"; 
		#print $#unique_male_nucleotides," ",$#unique_female_nucleotides,"\n";
		# looks fine
		if(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides != -1)){ #execute this block if the index number of the last element in the array is not -1 (which would only be true if the array were empty)
			# we can compare homoz and het genotypes because both sexes have data
			if(($#unique_male_nucleotides == 0)&&($#unique_female_nucleotides > 0)){ #if the number of the last index in @males is 0, then all males are homozyous at his position. if the number of the last index in @females is greater than 0, then at least some females are heretozygous
				# all males are homoz but at least some females are hets or homoz for another SNP
				# check if the proportion of divergent positions in females is high enough
				$diverged=0; #set diverged to 0 initially
				for ($x = 0 ; $x <= $#females ; $x++ ) { #x equals zero, x is less than or equal to the last index number of the females array, increment x by one each iteration
					if($females[$x] ne $males[0]){ #if the current female does not equal the first male (we can just use the first male, because in this case, we already know all males to be homozygous)
						$diverged+=1; #add 1 to the value of diverged. this effectively counts how many times a given divergent allele is found in the population of females.
					}
				}
				if($diverged > $proportion*($#females+1)){	#if the number of times an allele is divergent in the female population is greater than the proportion of divergence threshhold times the number of females, then the information for the position gets printed to the outfile: chrom 	position	"sex specific snp"	number of females	number of males
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_SNP\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_SNP
					# Category of 1 means ZW or a female-specific SNP
				}
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0; #start diverged at 0
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) { #the value of x is always in multiples of 2 so that you get both alleles for one individual
					if($females[$x] ne $females[$x+1]){ #if any female is heterozygous, add 1 to the value of diverged
						$diverged+=1;
					}
				}
				if($diverged == ($#females+1)/2){	#if the diverged count is equal to the number of females, print sex-specific positions to the outfile
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_heterozygosity\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_heterozygosity
					# 1 means all females are hets and all males are homoz
				}					
			}
			elsif(($#unique_male_nucleotides > 0)&&($#unique_female_nucleotides == 0)){
				# all females are homoz but at least some males are hets or homoz for another SNP
				# check if the proportion of divergent positions in males is high enough
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x++ ) {
					if($males[$x] ne $females[0]){
						$diverged+=1;
					}
				}
				if($diverged > $proportion*($#males+1)){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_SNP\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_SNP
					# -1 means XY or a male-specific SNP
				}	
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x=$x+2 ) {
					if($males[$x] ne $males[$x+1]){
						$diverged+=1;
					}
				}
				if($diverged == ($#males+1)/2){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_heterozygosity\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_heterozygosity
					# -1 means all males are hets and all females are homoz
				}
			}
			elsif(($#unique_male_nucleotides == 0)&&
			($#unique_female_nucleotides == 0)&&
			($unique_male_nucleotides[0] ne $unique_female_nucleotides[0])){
				# males are homoz, females are homoz, but fixed for a different diverged nucleotide
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tFixed_divergence\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
				# Fixed_divergence
				# 1 means diverged
			}
		}
		elsif(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides == -1)){ #male array is not empty, female array IS empty
			# females have no data
			# could be male-specific
			if((($#males +1)/2) > $proportion*$number_of_male_individuals_genotyped){
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_nucleotides\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n";
				# Sex_specific_nucleotides
				# -1 means male specific or male specific
			}	
		}
		elsif(($#unique_male_nucleotides == -1)&&($#unique_female_nucleotides != -1)){ #male array is empty, female array is not empty
			# males have no data
			# could be female-specific
			if((($#females +1)/2) > $proportion*$number_of_female_individuals_genotyped){
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_nucleotides\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
				# Sex_specific_nucleotides
				# 1 means fem specific or female specific
			}	
		}
	}
} # end while	
close OUTFILE1;

```
