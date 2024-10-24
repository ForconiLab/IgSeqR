# IgSeqR

IgSeqR is a tool designed for the identification, assembly, and characterization of full-length tumor immunoglobulin transcripts from unselected RNA sequencing data. This tool provides a user-friendly protocol for extracting and accurately determining the full tumor IG transcript sequence, which can be crucial for refining diagnosis, prognosis, and potential therapeutic targeting in B-cell tumors.

## Features

- Identification and assembly of full-length immunoglobulin transcripts.
- Characterization of tumor immunoglobulin sequences.
- Supports both BAM and FASTQ input formats.

## Installation

Clone the repository and run the setup script:
```bash
git clone https://github.com/ForconiLab/IgSeqR.git
cd IgSeqR
./setup.sh
```

### Key Tools and Versions

The following key tools and their versions are required for IgSeqR to function correctly:

- `blast=2.13.0`
- `hisat2=2.2.1`
- `kallisto=0.48.0`
- `samtools=1.16.1`
- `trinity=2.13.2`

The setup script will handle the installation of the Conda environment and necessary dependencies.

**Note:** This tool has been designed for and tested on Linux-based operating systems. Users on non-Linux systems can manually configure a conda environment using the provided dependencies, but this is beyond the scope of the tool.

## Usage

### Command Line Arguments

The script accepts the following command line arguments:

**Required:**

- `--bam, -b` : BAM file
- `--fastq, -fq` : Two FASTQ files (.fq, .fastq, .fq.gz, .fq.gz)
- `--hisat_ref, -hr` : HISAT reference

**Optional:**

- `--cores, -@` : Number of cores (default: number of available cores)
- `--mem, -m` : Memory in GB (default: total system memory)
- `--igh_ref` : IGH BLAST reference (.fasta, .fa)
- `--igkl_ref` : IGKL BLAST reference (.fasta, .fa)
- `--out, -o` : Output directory (default: ./igseqr)
- `--sample, -s` : Sample name (default: igseqr)
- `--chain, -c` : Chain type Heavy [H], Light [L], Both [B] (default)
- `--help` : Display help message

### Example Command

**BAM input:**
```bash
igseqr --bam sample.bam --hisat_ref hisat_index --out results --sample my_sample
```

**FASTQ input:**
```bash
igseqr --fastq sample_read1.fastq sample_read2.fastq --hisat_ref hisat_index --out results --sample my_sample
```

## References

For more detailed information about IgSeqR, please refer to the  [bioRxiv pre-print](https://www.biorxiv.org/content/10.1101/2024.09.03.611002v1).

## License

This project is licensed under the Apache License 2.0 - see the LICENSE file for details.

