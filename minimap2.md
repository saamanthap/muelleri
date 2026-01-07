Extract the first 500 entries from the Trinity assembly (so I can test on a smaller subset):
awk '/^>/ {n++} n>500 {exit} {print}' input.fasta > output.fasta
