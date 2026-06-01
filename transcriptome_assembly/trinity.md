I ran Trinity on Brian's machine. (Specifically info113, since it wasn't too busy at the time.) My script is called "trinity_script.sh":
```
#!/bin/bash

lin=/2/scratch/samp/newly_trimmed_reads/*R1.fq.gz
rin=/2/scratch/samp/newly_trimmed_reads/*R2.fq.gz

ll=$(
        for l in ${lin}; do
                basename $l
        done
)

rr=$(
        for r in ${rin}; do
                basename $r
        done
)

lstring=$(echo ${ll} | sed 's/ /,/g')
rstring=$(echo ${rr} | sed 's/ /,/g')

/opt/local/trinity/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left ${lstring} --right ${rstring} --max_memory 50G --CPU 6 --output sam_trinity_assembly

```
I used this command, from the directory that contains the trimmed reads ("newly_trimmed_reads"):
```
 ../trinity_script.sh 1> trinity.out 2> trinity.err &
```

