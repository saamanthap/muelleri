In order to check if primers will bind to a SINGLE site in the genome (with respect to the borealis reference), you need to blast the primer sequences against the borealis database on commandline. Here is an example command: 
```
blastn -query rbbp4_primers.fa -db ../../Austin_genome_Xborealis/Xbo.v1.fa_blastable -outfmt 6 -task blastn-short > rbbp4_primers.out
```
Now take the output file and check that there is only ONE hit that covers the full query length and has 100% percent identity. You can use sort to glance at this for each primer. (I used grep to easily search for the primer name.)
```
sort -n -r -s -k4,2 -k3,1 rbbp4_primers.out | grep '1_F' | head
sort -n -r -s -k4,2 -k3,1 rbbp4_primers.out | grep '1_R' | head
```
Example output: 
```
1_R     Chr2S   100.000 20      0       0       1       20      8054298980543008        0.007   40.1
1_R     Chr2S   100.000 16      0       0       1       16      9256657492566559        1.7     32.2
1_R     Chr2S   94.737  19      1       0       2       20      7463259374632575        6.6     30.2
```
In the above, there are multiple primers with 100% percent identity, but only one that covers the full primer length. This is good! You can fiddle with the annealing temp to ensure specificity when using the primers.
