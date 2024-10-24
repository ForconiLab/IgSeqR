#!/bin/bash

# This script is designed to download and extract the HISAT indexed reference transcriptome
# of the GRCh38 Homo sapiens genome build. The reference transcriptome is provided as a 
# compressed tarball (.tar.gz) file, which combines both archiving and compression.
# 
# The script performs the following steps:
# 1. Creates the output directory if it does not already exist.
# 2. Downloads the compressed tarball from the specified URL.
# 3. Verifies the success of the download.
# 4. Extracts the contents of the tarball into the output directory.
# 5. Verifies the success of the extraction.
#
# Usage:
# ./get_hisat_indexed_ref.sh
#
# Ensure that you have the necessary permissions to create directories and write files
# in the specified output location. Additionally, make sure that wget and tar are installed
# on your system.

# Define the URL and output directory
URL="https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz"
OUTPUT_FILE="grch38_snptran.tar.gz"

# Download the file
echo -e "Downloading $URL...\n"
wget $URL -O $OUTPUT_FILE

# Check if the download was successful
if [ $? -ne 0 ]; then
    echo " - ERROR: Download failed!"
    exit 1
fi

# Extract the tar.gz file
echo "Extracting $OUTPUT_FILE..."
tar -xzvf $OUTPUT_FILE

# Check if the extraction was successful
if [ $? -ne 0 ]; then
    echo " - ERROR: Extraction failed!"
    exit 1
fi

echo "Download and extraction completed successfully."
