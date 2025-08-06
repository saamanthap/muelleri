Get the acronym after GenePageIDs: 
```
awk '{print gensub(/.*GenePageIDs:([^|]+)\|.*/, "\\1", "g")}' m_f_laevis_transcriptome_blast_results.out > m_f_annotations
```
Find a string inside () and replace the entire line with the string. This is a command for inside vi: 
```
%s/.*(\(.*\)).*/\1/
```
