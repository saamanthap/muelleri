# Subsetting muel_res.txt
muel_res.txt contains the output of the DESeq2 analysis of differential expression. The only information that has been pre-filtered is transcripts with a count sum of less than 2 for all samples.  

The file looks like this: 
```
"baseMean"      "log2FoldChange"        "lfcSE" "stat"  "pvalue"        "padj"
"TRINITY_DN37947_c0_g1_i1"      2.32589035223894        -1.10133580543906       1.0727632257983 -1.02663456292463       0.304592572718961       NA
"TRINITY_DN816_c0_g1_i2"        19.1745585951407        0.102368878209317       0.694342719096557       0.147432781238801       0.8827904325416 0.999986614612765
"TRINITY_DN17376_c0_g1_i14"     17.5829886832452        -1.60326581098384       1.74816592907749        -0.917113063649452      0.359083396197561       0.999986614612765
"TRINITY_DN14109_c2_g2_i1"      2.1227073860638 1.271017635168  1.28469026599311        0.989357255062146       0.322488380296208       NA
"TRINITY_DN23737_c1_g4_i1"      6.86482938282486        -0.107900175060747      0.605669940780386       -0.178150124012628      0.858605078312968       0.999986614612765
"TRINITY_DN9013_c1_g1_i21"      1846.56982351543        0.263657099599296       0.450085799901483       0.585792974710614       0.558014645330457       0.999986614612765
"TRINITY_DN7344_c0_g1_i3"       4.75232726704878        -5.14430409474999       3.36740068595493        -1.52767804443538       0.126592493930208       0.999986614612765
"TRINITY_DN11162_c0_g2_i3"      4.10513861913165        1.07843868312973        1.6348334569038 0.659662718899929       0.509470296024432       0.999986614612765
"TRINITY_DN6591_c0_g1_i9"       10.5436807514657        -1.09354687729788       1.86950548023553        -0.584939112968051      0.558588657852492       0.999986614612765
```
Fixing the file so that all values are separated by tabs:
```
awk -v OFS='\t' '{$1=$1; print}' your_file.txt > fixed_file.txt
```




Note that not all lines are perfectly tab-delimited (not sure why). Most of the time, I want to filter by log2FoldChange (the third column) and pvalue (the sixth column). The headings row is treated as a regular row, so be sure to exclude or ignore it.  

Negative log2FoldChange values imply female-biased expression. Positive log2FoldChange values imply male-biased expression.
