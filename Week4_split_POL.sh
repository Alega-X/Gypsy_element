#!/bin/bash

# Use a for loop to go through all of POL gene sequences in week 3 results
for i in {1..1378}; do
    # Check whether the line number is an even number
    if [ $((i % 2)) -eq 0 ]; then
        # Extract and save each POL gene into an independent fasta file for HHpred job submission
        LINE=$(awk -v line_num=$(($i - 1)) 'NR==line_num {gsub(/>/, ""); print}' POL_result.fasta)
        awk "NR==$(($i - 1))" POL_result.fasta >> "/scratch/lf10/zx6715/HHdata/${LINE}.fasta"
        awk "NR==$i" POL_result.fasta >> "/scratch/lf10/zx6715/HHdata/${LINE}.fasta"
    fi
done
