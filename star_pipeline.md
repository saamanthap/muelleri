# Star Pipeline 
I'll use this to  
* Map X. muelleri RNA-Seq reads to the X. laevis reference genome
* Get counts for each gene
* Perform DESeq2 analysis of differential expression

The RNA-Seq reads are stored here: 
```
/home/samp/projects/rrg-ben/for_Sam/muel/muel_fq/
```
The file names look like this, although there are other files in this folder: 
```
X_muelleri_tad31_S11_L001__trim_R1.fq.gz
X_muelleri_tad31_S11_L001__trim_R2.fq.gz
X_muelleri_tad32_S12_L001__trim_R1.fq.gz
X_muelleri_tad32_S12_L001__trim_R2.fq.gz
```
I'm using the X. laevis reference genome 10.1 (which can be downloaded from Xenbase using wget command). I need the associated annotation file in .gtf format: 
```
/home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/XENLA_10.1_Xenbase.gtf
```
## Star Mapping 
Path to the script: 
```
/home/samp/projects/rrg-ben/for_Sam/muel/ben_scripts/sam_2025_STAR_mapreads.sh
```
The script: 
(Note that it is important to include the --sjdbGTFfile option because it should improve the accuracy of mapping - I was having a problem before where 50% of reads are unmapped because they are too short.)
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=64gb
#SBATCH --output=STAR_map.%J.out
#SBATCH --error=STAR_map.%J.err
#SBATCH --account=rrg-ben
#SBATCH --mail-user=pottss5@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load StdEnv/2023 star/2.7.11b

STAR --genomeDir /home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/ \
--sjdbGTFfile /home/samp/projects/rrg-ben/for_Sam/2021_XL_v10_refgenome/XENLA_10.1_Xenbase.gtf \
--runThreadN 6 \
--readFilesIn ${1} ${2} \
--outFileNamePrefix ${3} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--alignEndsProtrude 10 ConcordantPair \
        --readFilesCommand zcat
```
Make sure to check the outfile ending with "Log.final.out". An example of one of my jobs that mapped only 50% of reads (not ideal): 
```
  Started job on |       Jun 13 22:01:06
                             Started mapping on |       Jun 13 22:01:28
                                    Finished on |       Jun 14 00:49:33
       Mapping speed, Million of reads per hour |       8.44

                          Number of input reads |       23656227
                      Average input read length |       279
                                    UNIQUE READS:
                   Uniquely mapped reads number |       9512813
                        Uniquely mapped reads % |       40.21%
                          Average mapped length |       259.74
                       Number of splices: Total |       13713870
            Number of splices: Annotated (sjdb) |       13362769
                       Number of splices: GT/AG |       13587623
                       Number of splices: GC/AG |       72983
                       Number of splices: AT/AC |       6728
               Number of splices: Non-canonical |       46536
                      Mismatch rate per base, % |       3.16%
                         Deletion rate per base |       0.06%
                        Deletion average length |       2.73
                        Insertion rate per base |       0.05%
                       Insertion average length |       2.74
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       953597
             % of reads mapped to multiple loci |       4.03%
        Number of reads mapped to too many loci |       112194
             % of reads mapped to too many loci |       0.47%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |       0
       % of reads unmapped: too many mismatches |       0.00%
            Number of reads unmapped: too short |       13046141
                 % of reads unmapped: too short |       55.15%
                Number of reads unmapped: other |       31482
                     % of reads unmapped: other |       0.13%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
```

