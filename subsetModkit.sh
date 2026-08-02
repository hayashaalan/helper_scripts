#!/bin/bash

# This is a script to filter and write modkit read-level data
# I created this script to make the data compatible with NanoMethViz; running the entire file through the tool is crazy
# Provide the chromosome throught the commandline e.g src/modkit_filter.sh chr18
chrom="$1"
input_dir="results/gambl/modkit-1.0/99-outputs/tsv/promethION--hg38/"
output_dir="results/gambl/modkit-1.0/level3/${chrom}/"

mkdir -p "$output_dir"

for file in "$input_dir"/*.tsv.gz; do
    filename=$(basename "$file")
    output_file="$output_dir/$filename"

    # Pull header and any lines with chr passed through commandline
    # Filter for mod_qual = m
    zcat "$file" | awk -v chrom="$chrom" '/^chrom/ {print; found=1; next} found && $1 == chrom && $4 == "m"' | gzip > "$output_file"
done

