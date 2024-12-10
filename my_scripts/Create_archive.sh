#!/bin/bash

# Define the target directory
TARGET_DIR="AMR_archive"

# Create the target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Define the name of the file to be moved
FILE_NAME="bwt_out/bwt.artifacts_mapping_stats.txt"

# Find all subdirectories and move the specified file
#find samples_2/ -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
#    if [ -f "$dir/$FILE_NAME" ]; then
        #echo "$dir/$FILE_NAME"
        # Get the base name of the directory
        #base_name=$(basename "$dir")
        #echo "$base_name"
        # Move and rename the file
        #echo "Moving $dir/$FILE_NAME .....to..... $TARGET_DIR/${base_name}_${FILE_NAME}"  
        #mv "$dir/$FILE_NAME" "$TARGET_DIR/${base_name}_${FILE_NAME}"
#    fi
#done

#############################################################################

TARGET_DIR="AMR_archive"

# Create the target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Define the name of the file to be moved
FILE_NAME="reads_2.fastq"

# Find all subdirectories and move the specified file
find samples_2/ -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
    if [ -f "$dir/$FILE_NAME" ]; then

        base_name=$(basename "$dir")

        echo "Moving $dir/$FILE_NAME .....to..... $TARGET_DIR/${base_name}_R2.fastq"
        #echo "Moving $dir/$FILE_NAME .....to..... $TARGET_DIR/${base_name}.fastq"
        mv "$dir/$FILE_NAME" "$TARGET_DIR/${base_name}_R2.fastq"
    fi
done