#!/usr/bin/env perl
use strict;
use warnings;

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
my @columns;
my @pat;
my @pat1;
my $x;

my @sexes = split("",$sexes);

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2  $!\n\n";
	exit;
}
print "Creating output file: $outputfile2\n";



if ($vcf =~ /.gz$/) {
	#open DATAINPUT, '<:gzip', $vcf or die "Could not read from $vcf: $!";
	open(DATAINPUT, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
}
else {
	open DATAINPUT, $vcf or die "Could not read from $vcf: $!";
}

my $number_of_samples=0;
my $switch1=0;
my $switch2=0;
my $num_daughters=0;
my $num_het_daughters=0;
my $num_sons=0;
my $num_het_sons=0;


while ( my $line = <DATAINPUT>) {
	@columns=split("	",$line);
		#print $line,"\n";
		if($columns[0] =~ m/^#/){ # this is a commented line
			if($columns[0] eq '#CHROM'){ # this is the first line
				$number_of_samples = scalar(@columns)-9;
				print "Number of samples ",$number_of_samples,"\n";
			}
		}
		else{ # this is the genotype data
			$switch1=0;
			$switch2=0;
			$num_daughters=0;
			$num_het_daughters=0;
			$num_sons=0;
			$num_het_sons=0;
			for ($x = 1 ; $x <= $number_of_samples; $x++ ){ # cycle through each sample
				if($sexes[$x-1] != 2){ # only consider samples that are included
					@pat = split(":",$columns[$x+8]); # this is the whole genotype and other info from the vcf file
					@pat1 = split(/[\|\/]/,$pat[0]); # this is only the genotype
					# check if the daughters are homozygous or missing
					if(($pat[0] ne './.')&&($pat[0] ne '.|.')&&
					($pat[0] ne '.')&&($sexes[$x-1] == 1)){ # this is a daughter with a genotype
					# this is a daughter
					$num_daughters+=1;
						if(($pat1[0] eq $pat1[1])&&($switch1 != 9)){
							# this daughter is homoz 
							$switch1 = 1;
							#print "goodpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
						}
						elsif(($pat1[0] ne $pat1[1])&&($switch1 != 9)){
							# this daughter is the first heteroz
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$switch1 = 9;
							$num_het_daughters+=1;
						}	
						elsif(($pat1[0] ne $pat1[1])&&($switch1 == 9)){
							# this daughter is heteroz, but not the first one
							#print "badpatd ",$columns[0],"\t",$columns[1]," ",$pat[0]," ",$pat[1]," ";
							$num_het_daughters+=1;
						}	
					}
					# check if the sons are heterozygous or missing
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
						elsif(($pat1[0] ne $pat1[1])&&($switch1 == 9)){
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
			}
			
		} # end else
} # end while
close DATAINPUT;
close OUTFILE;
close OUTFILE2;
