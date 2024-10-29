# Tutorial: Accessing FASTQ Reads from ArrayExpress/BioStudies for IgSeqR Pipeline

This tutorial will guide you through the steps to access FASTQ reads from the ArrayExpress/BioStudies database for the Chronic Lymphocytic Leukemia (CLL) sample with accession number E-MTAB-12017. We will use sample 822 as an example dataset for the IgSeqR pipeline.

## Overview

In this tutorial, you will learn how to:
1. Set up your environment for running the IgSeqR pipeline.
2. Download necessary data files from ArrayExpress/BioStudies.
3. Run the IgSeqR pipeline to analyse immunoglobulin sequences.

## Prerequisites

Before you begin, ensure you have access to a Linux-based terminal or command line interface. This is essential for running the commands and scripts provided in this tutorial. If you are using a different operating system, consider using a virtual machine or a Docker container with a Linux environment.

**Note:** The setup and use of virtual machines or Docker containers are beyond the scope of this tutorial. For more information, you can refer to the following resources:

- [Virtual Machines](https://www.freecodecamp.org/news/what-is-a-virtual-machine-and-how-to-setup-a-vm-on-windows-linux-and-mac/)
- [Docker](https://docker-curriculum.com/)

You will also need:
- An internet connection
- `wget` or `curl`, `git` and `conda` installed on your system

### Setting Up Your Environment

Clone the IgSeqR github repository to your environment and run the setup script:

```bash
git clone https://github.com/ForconiLab/IgSeqR.git
cd IgSeqR
./setup.sh
```

The following key tools and their versions are required for IgSeqR to function correctly:

- `blast=2.13.0`
- `hisat2=2.2.1`
- `kallisto=0.48.0`
- `samtools=1.16.1`
- `trinity=2.13.2`

The setup script will handle the installation of the Conda environment and necessary dependencies.

## Step 1: Create a Directory for Analysis

First, create and navigate to a new directory to conduct the analysis of sample 822:

```bash
mkdir 822_igseqr
cd 822_igseqr
```

## Step 2: Access ArrayExpress/BioStudies

1. Open your web browser and navigate to the ArrayExpress website.
2. In the search bar, enter the accession number `E-MTAB-12017` and press Enter.
3. Click on the study title to open the study page.

## Step 3: Download FASTQ Files

On the study page, you will find links to the raw data files. Follow these steps to download the FASTQ files for sample 822:

1. Locate the table with the raw data files under the heading "Detailed sample information and links to data".
2. Find the sample(s) of interest and locate the link to the ENA submission.
3. Click the download link for the submitted FASTQ file.
4. Save these files to the newly created sample directory.

Alternatively, you can use the following commands in your terminal to download the files for sample 822 directly:

### Using `wget`

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10015005/WTCHG_433964_261_1_trimmed.fastq.gz -O 822_read1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10015006/WTCHG_433965_261_1_trimmed.fastq.gz -O 822_read2.fastq.gz
```

### Using `curl`

```bash
curl -o 822_read1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10015005/WTCHG_433964_261_1_trimmed.fastq.gz
curl -o 822_read2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR100/ERR10015006/WTCHG_433965_261_1_trimmed.fastq.gz
```

For more information, refer to the [ArrayExpress documentation](https://www.ebi.ac.uk/biostudies/arrayexpress/help).

## Step 4: Download Hisat Index

The Hisat indexed reference transcriptomes for various genome builds can be found [here](https://daehwankimlab.github.io/hisat2/download/#h-sapiens).

First, create a new folder to store the reference files:

```bash
mkdir -p HISAT_REF
```

Next, download the GRCh38 genome_snp_tran build:

### Using `wget`

```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz -O ./HISAT_REF/grch38_snptran.tar.gz
```

### Using `curl`

```bash
curl -o ./HISAT_REF/grch38_snptran.tar.gz https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz
```

Uncompress the downloaded archive:

```bash
tar -xzvf ./HISAT_REF/grch38_snptran.tar.gz -C ./HISAT_REF/
```

**Note:** This download is approximately >4GB and may take some time. Ensure you have sufficient disk space and a stable internet connection.

Alternatively, the `get_hisat_indexed_ref.sh` script located in the `data/HISAT` section of this repository can be run to perform this step.

## Step 5: Run IgSeqR

With the two FASTQ files and the Hisat indexed reference, you can now analyze the data using the IgSeqR pipeline:

```bash
conda activate igseqr

igseqr --fastq 822_read1.fastq.gz 822_read2.fastq.gz \
--hisat_ref ./HISAT_REF/grch38_snp_tran/genome_snp_tran \
--out igseqr \
--sample CLL822 \
--chain B
```

Here:
- `--fastq` supplies the two read FASTQ files.
- `--hisat_ref` specifies the location of the Hisat index files. Note that the format must be the path and basename of the indexed reference transcriptome files, if step 4 was followed this will be `./HISAT_REF/grch38_snp_tran/genome_snp_tran`.
- `--out` is the directory where the results will be stored.
- `--sample` is the sample name or identifier.
- `--chain` specifies the immunoglobulin chain(s) to analyze. `B` has been selected for both heavy and light chain analysis.

## Conclusion

Upon completion of this tutorial, you should have generated the five most abundant immunoglobulin transcripts for both the heavy and the light chain. The resultant FASTA files can be uploaded to the IMGT/V-QUEST webtool or other immunogenetics annotation tools for downstream analysis. For more detailed information about IgSeqR, please refer to the [bioRxiv pre-print](https://www.biorxiv.org/content/10.1101/2024.09.03.611002v1).