#!/bin/bash
# Use a for loop to go through all gypsy element CDS
for i in {1..2668}; do
    # Extract the i-th line from the FASTA file
    line=$(awk "NR==$i" Drosophila_gypsy_RepBase_2024-Aug.CDS.fasta)
    
    # Check if the line is present in POL_real.txt
    if grep -Fq "$line" POL_real.txt; then
        # Extract the i-th line from the FASTA file and write it to POL_result.fasta
        awk "NR==$i" Drosophila_gypsy_RepBase_2024-Aug.CDS.fasta >> POL_result.fasta
        
        # Extract the next line (i+1) from the FASTA file and append it to POL_result.fasta
        awk "NR==$(($i + 1))" Drosophila_gypsy_RepBase_2024-Aug.CDS.fasta >> POL_result.fasta
    fi
done
