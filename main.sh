#!/bin/bash

# Define the file paths
chr_length_file="chromelength.txt"
gene_gtf_file="Homo_sapiens.GRCh37.87.gtf"
input_file="rpkm.txt"

# Define step size
step_size=1

# Loop through width values with the specified step size
for ((width = 1; width <= 50; width += step_size))
do
    echo "Running with width = $width"
    python main.py --chr_length_file $chr_length_file --gene_gtf_file $gene_gtf_file --width $width --input_file $input_file
done
