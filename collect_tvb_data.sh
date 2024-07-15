#!/bin/bash

# Find and copy weights.csv and lengths.csv files
for pattern in "sub-*_V1_MR_parc-mmp1_sc_weights.csv" "sub-*_V1_MR_parc-mmp1_sc_lengths.csv"; do
    find /Volumes/bronkodata_work/hcpython/testsub -path "*/T1w/Diffusion/MRTrix/${pattern}" | while read file; do
        # Extract the subject ID from the file path
        subjectid=$(basename $file | cut -d'_' -f1 | cut -d'-' -f2)

        # Create the target directory if it doesn't exist
        target_dir="/Volumes/bronkodata_work/hcpython/testsub/SCs/${subjectid}"
        mkdir -p "$target_dir"

        # Copy the file to the target directory
        cp "$file" "$target_dir"
    done
done