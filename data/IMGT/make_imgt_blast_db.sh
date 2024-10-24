#!/bin/bash

# This script downloads reference files for the Heavy and Light chains of the Homo sapiens
# immunoglobulin (IG) genes from the IMGT database. It then processes these files to create
# combined FASTA files for each chain and generates BLAST databases from these combined files.
#
# The script performs the following steps:
# 1. Creates directories for the Heavy (IGH) and Light (IGKL) chains.
# 2. Downloads the relevant FASTA files for each chain from the IMGT database.
# 3. Combines the downloaded FASTA files into single files (IGH.fasta and IGKL.fasta) for each chain.
# 4. Processes the sequences to remove dots and wrap lines at 60 characters.
# 5. Creates BLAST databases from the combined FASTA files.
#
# Usage:
# ./make_imgt_blast_db.sh
#
# Ensure that you have the necessary permissions to create directories and write files
# in the specified output locations. Additionally, make sure that wget and makeblastdb
# are installed on your system.

set -e

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for required commands
for cmd in wget makeblastdb; do
    if ! command -v $cmd >/dev/null 2>&1; then
        echo "ERROR: $cmd is not installed, please install blast or activate the igseqr conda environment."
        exit 1
    fi
done

# Function to download a file and check for errors
download_file() {
    local url=$1
    local output=$2
    echo "Downloading $url..."
    wget $url -O $output
    if [ $? -ne 0 ]; then
        echo "Failed to download $url"
        exit 1
    fi
}

# Function to process FASTA files
process_fasta() {
    local input_files=("$@")
    local output_file=${input_files[-1]}
    unset input_files[-1]

    echo "" > $output_file

    cat "${input_files[@]}" | while read line; do
        if [[ $line == \>* ]]; then
            # write previous sequence and header to file
            if [[ -n $sequence ]]; then
                echo ">$header" >> $output_file
                # remove dots from sequence and wrap lines at 60 characters
                echo "$sequence" | tr -d '.' | fold -w 60 >> $output_file
                sequence=""
            fi
            header=$(echo $line | cut -d "|" -f 2)
            sequence=""
        else
            sequence+="$line"
        fi
    done

    # write last sequence and header to file if not empty
    if [[ -n $sequence ]]; then
        echo ">$header" >> $output_file
        echo "$sequence" | tr -d '.' | fold -w 60 >> $output_file
    fi
}


# Function to create BLAST database
create_blast_db() {
    local fasta_file=$1
    local db_title=$2
    makeblastdb -in $fasta_file \
    -input_type fasta \
    -parse_seqids \
    -title "$db_title" \
    -dbtype nucl
}

## Heavy Chain
mkdir -p IGH
cd IGH

download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta" "IGHV.fasta"
download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta" "IGHD.fasta"
download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta" "IGHJ.fasta"

process_fasta IGHV.fasta IGHD.fasta IGHJ.fasta IGH.fasta
create_blast_db IGH.fasta "IGH_db"

cd ..

## Light Chain
mkdir -p IGKL
cd IGKL

download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta" "IGKV.fasta"
download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKJ.fasta" "IGKJ.fasta"
download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta" "IGLV.fasta"
download_file "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLJ.fasta" "IGLJ.fasta"

process_fasta IGKV.fasta IGKJ.fasta IGLV.fasta IGLJ.fasta IGKL.fasta
create_blast_db IGKL.fasta "IGKL_db"

cd ..

echo "All tasks completed successfully."
