#!/bin/bash

################################## Parse Arguments for Configuration Variables ##################################

# Initialise variables with default values
CORES=$(nproc)
TOTAL_SYS_MEM=$(awk '/MemTotal/ {printf "%.0f\n", $2/1024/1024}' /proc/meminfo)
MEM=$TOTAL_SYS_MEM
OUT_DIR="./igseqr"
SAMPLE="igseqr"
CHAIN_TYPE="B"
FASTQ=()
HOME_DIR=$(pwd)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KEEP_INT=false

# Function to display help message
show_help() {
    echo "Usage: $0 [options]"
    echo "Required:"
    echo "  --bam, -b          BAM file"
    echo "  --fastq, -fq       Two FASTQ files"
    echo "  --hisat_ref, -hr   HISAT reference"
    echo "Optional:"
    echo "  --cores, -@        Number of cores (default: $(nproc))"
    echo "  --mem, -m          Memory in GB (default: $TOTAL_SYS_MEM)"
    echo "  --igh_ref          IGH reference"
    echo "  --igkl_ref         IGKL reference"
    echo "  --out, -o          Output directory (default: ./igseqr)"
    echo "  --sample, -s       Sample name (default: igseqr)"
    echo "  --chain, -c        Chain type Heavy [H], Light [L], Both [B] (default)"
    echo "  --keep_int         Keep all intermediate files"
    echo "  --help             Display this help message"
}

# Function to check if a file exists and is of the correct format
check_file() {
    local file=$1
    local format=$2
    if [[ ! -e "$file" ]]; then
        echo "ERROR: The specified file '$file' does not exist."
        exit 1
    elif [[ ! "$file" =~ \.$format(\.gz)?$ ]]; then
        echo "ERROR: The specified file '$file' should be in '$format' format."
        exit 1
    fi
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam|-b) BAM="$2"; shift 2 ;;
        --fastq|-fq)
            if [[ -n "$2" && -n "$3" && "$2" != -* && "$3" != -* ]]; then
                FASTQ+=("$2" "$3")
                shift 3
            else
                echo "ERROR: --fastq requires exactly two FASTQ files."
                exit 1
            fi
            ;;
        --cores|-@) CORES="$2"; shift 2 ;;
        --mem|-m) MEM="$2"; shift 2 ;;
        --hisat_ref|-hr) HISAT_REF="$2"; shift 2 ;;
        --igh_ref) IGH_REF="$2"; shift 2 ;;
        --igkl_ref) IGKL_REF="$2"; shift 2 ;;
        --out|-o) OUT_DIR="$2"; shift 2 ;;
        --sample|-s) SAMPLE="$2"; shift 2 ;;
        --chain|-c) CHAIN_TYPE=$(echo "$2" | tr '[:lower:]' '[:upper:]' | cut -c1); shift 2 ;;
        --keep_int) KEEP_INT=true; shift ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
done

# Check for required arguments
if [[ -z "$BAM" && ${#FASTQ[@]} -eq 0 ]]; then
    echo "ERROR: Either --bam or --fastq must be provided."
    show_help
    exit 1
elif [[ -n "$BAM" && ${#FASTQ[@]} -gt 0 ]]; then
    echo "ERROR: Both BAM and FASTQ files are provided. Please provide only one type."
    show_help
    exit 1
elif [[ -n "$BAM" ]]; then
    check_file "$BAM" "bam"
else
    for fq in "${FASTQ[@]}"; do
        check_file "$fq" "fastq|fq"
    done
    # Assign FASTQ files to read1 and read2
    READ1="${FASTQ[0]}"  # First FASTQ file
    READ2="${FASTQ[1]}"  # Second FASTQ file
fi

# Confirm HISAT reference has been provided
if [[ -z "$HISAT_REF" ]]; then
    echo "ERROR: --hisat_ref is required."
    show_help
    exit 1
fi

# Check if MEM is an integer
if ! [[ "$MEM" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Memory argument must be an integer."
    exit 1
elif (( MEM > TOTAL_SYS_MEM )); then
    echo "ERROR: Memory argument (${MEM}) exceeds total system memory (${TOTAL_SYS_MEM} GB)."
    exit 1
fi

# Validate CHAIN_TYPE
case $CHAIN_TYPE in
    H) CHAIN_TYPE="HEAVY" ;;
    L) CHAIN_TYPE="LIGHT" ;;
    B) CHAIN_TYPE="BOTH" ;;
    *) echo "ERROR: Invalid chain type. Specify Heavy [H], Light [L], or Both [B] (default)."; exit 1 ;;
esac

if [[ "$CHAIN_TYPE" == "HEAVY" || "$CHAIN_TYPE" == "BOTH" ]]; then
    IG_CHAIN+=(IGH)
    if [[ -n "$IGH_REF" ]]; then
        check_file "$IGH_REF" "fasta|fa"
    else
        IGH_REF="$SCRIPT_DIR/data/igseqr/IMGT/IGH/IGH.fasta"
    fi
fi

if [[ "$CHAIN_TYPE" == "LIGHT" || "$CHAIN_TYPE" == "BOTH" ]]; then
    IG_CHAIN+=(IGKL)
    if [[ -n "$IGKL_REF" ]]; then
        check_file "$IGKL_REF" "fasta|fa"
    else
        IGKL_REF="$SCRIPT_DIR/data/igseqr/IMGT/IGKL/IGKL.fasta"
    fi
fi

mkdir $OUT_DIR

# If BAM file is provided
if [[ -n $BAM ]]; then
    echo 
    echo ===============================================================================
    echo IgSeqR - EXTRACT TRANSCRIPTS FROM BAM
    echo ===============================================================================
    echo
    echo ================================= Config ======================================
    echo "BAM:                    $BAM"
    echo "Sample ID:              $SAMPLE"
    echo "Hisat Reference:        $HISAT_REF"
    echo "IG Chain(s):            ${IG_CHAIN[@]}"
    for element in "${IG_CHAIN[@]}"; do
        if [[ "$element" == "IGH" ]]; then
            echo "IG Heavy Reference:     $IGH_REF"
        elif [[ "$element" == "IGKL" ]]; then
            echo "IG Light Reference:     $IGKL_REF"
        fi
    done
    echo "Cores:                  $CORES"
    echo "Memory:                 $MEM"
    echo "Output Directory:       $OUT_DIR"
    echo =============================================================================== 

    #################################### Convert Bam into FASTQ #################################### 
    ## Convert BAM reads back into FASTQ for HISAT2 alignment
    echo
    date
    echo =========================== Converting BAM to FASTQ ===========================
    FASTQ_DIR=$OUT_DIR/FASTQ
    mkdir $FASTQ_DIR

    echo Sorting BAM by name ... 
    samtools sort -n -@$CORES $BAM -o $FASTQ_DIR/"$SAMPLE"_sorted.bam

    echo Converting BAM to FASTQ ...
    samtools fastq -@ $CORES -n -c 6 $FASTQ_DIR/"$SAMPLE"_sorted.bam \
    -1 $FASTQ_DIR/"$SAMPLE"_read1.fastq.gz \
    -2 $FASTQ_DIR/"$SAMPLE"_read2.fastq.gz \
    -0 /dev/null -s /dev/null

    FASTQ1="$FASTQ_DIR/"$SAMPLE"_read1.fastq.gz"
    FASTQ2="$FASTQ_DIR/"$SAMPLE"_read2.fastq.gz"
fi

# If FASTQ files are provided, check if they exist
if [[ -n $READ1 && -n $READ2 ]]; then
    echo 
    echo ===============================================================================
    echo IgSeqR - EXTRACT TRANSCRIPTS FROM FASTQ
    echo ===============================================================================
    echo
    echo ================================= Config ======================================
    echo "FASTQ Read 1:           $READ1"
    echo "FASTQ Read 2:           $READ2"
    echo "Sample ID:              $SAMPLE"
    echo "Hisat Reference:        $HISAT_REF"
    echo "IG Chain(s):            ${IG_CHAIN[@]}"
    for element in "${IG_CHAIN[@]}"; do
        if [[ "$element" == "IGH" ]]; then
            echo "IG Heavy Reference:     $IGH_REF"
        elif [[ "$element" == "IGKL" ]]; then
            echo "IG Light Reference:     $IGKL_REF"
        fi
    done
    echo "Cores:                  $CORES"
    echo "Memory:                 $MEM"
    echo "Output Directory:       $OUT_DIR"
    echo ===============================================================================

    FASTQ1=$READ1
    FASTQ2=$READ2
fi

######################################## HISAT2 Alignment ########################################
## Map to genome using hisat2, sort and index bam
HISAT_DIR=$OUT_DIR/HISAT
mkdir $HISAT_DIR

echo
date
echo ========================= Performing HISAT2 Alignment =========================

hisat2 -p $CORES --phred33 -t -x $HISAT_REF -1 $FASTQ1 -2 $FASTQ2 | \
  samtools view -@$CORES -bS -o - - | \
  samtools sort -@$CORES - -o $HISAT_DIR/hisat2_output.bam &&
samtools index -@$CORES $HISAT_DIR/hisat2_output.bam -o $HISAT_DIR/hisat2_output.bam.bai

echo
echo HISAT2 Flagstat:
samtools flagstat $HISAT_DIR/hisat2_output.bam

############################# Clean FASTQ Files if created  #####################################
if ! $KEEP_INT; then
    if [[ -e "$FASTQ_DIR" ]]; then
        rm -r $FASTQ_DIR
    fi
fi

######################################### Read Selection ######################################### 
## Select unmapped reads and any reads mapping to IG loci to take forward to trinity
## Remove anything not required to assemble IG transcripts
## The fewer reads, the better as trinity is very slow

## Coordinates to keep with buffer:
##  - IGH: "14:100000000-110000000"
##  - IGK: "2:87000000-92000000"
##  - IGL: "22:20500000-24500000"

## Coordinates to keep without buffer:
##  - IGH: "14:105550000-106900000" (IGHA2 to IGHVIII-82)
##  - IGK: "2:88697000-92240000" (IGKC to IGKV1OR2 pseudogene)
##  - IGL: "22:22005000-23590000" (IGLV4-69 to IGLL1)

echo
date
echo ============================== IG Read Filtering ==============================

IGH="14:105550000-106900000"
IGK="2:88697000-92240000"
IGL="22:22005000-23590000"
REGIONS=""

echo "Filtering for unmapped reads and reads mapping to:" 

# Construct the regions to extract based on CHAIN
for element in "${IG_CHAIN[@]}"; do
    if [[ "$element" == "IGH" ]]; then
        REGIONS="${REGIONS}$IGH "
        echo " - IGH: $IGH"
    elif [[ "$element" == "IGKL" ]]; then
        REGIONS="${REGIONS}$IGK $IGL "
        echo " - IGK: $IGK"
        echo " - IGL: $IGL"
    fi
done

# extract unmapped reads and reads mapping to specified regions and sort the bam by name for fastq extraction
samtools merge -f - \
    <(samtools view -@ $CORES -b -f 4 $HISAT_DIR/hisat2_output.bam) \
    <(samtools view -@ $CORES -b $HISAT_DIR/hisat2_output.bam $REGIONS) | \
samtools sort -n -@$CORES -o $HISAT_DIR/keep_reads.bam


#################################### Convert Back into FASTQ ################################### 
## Convert filtered reads back into FASTQ for running through Trinity

echo
date
echo =========================== Converting BAM to FASTQ ===========================

IG_READS_DIR=$OUT_DIR/IG_READS
mkdir $IG_READS_DIR

samtools fastq -@ $CORES -n -c 6 $HISAT_DIR/keep_reads.bam \
-1 $IG_READS_DIR/"$SAMPLE"_read1.fastq.gz \
-2 $IG_READS_DIR/"$SAMPLE"_read2.fastq.gz \
-0 /dev/null -s /dev/null

###################################### Clean HISAT2 Files #####################################
if ! $KEEP_INT; then
    if [[ -e "$HISAT_DIR" ]]; then
        rm -r $HISAT_DIR
    fi
fi

########################### Trinity genome naive trascript assembly ###########################
## here we assemble the kept reads into transcripts. 
## Option choice is important, please note the following:
##   - no digital normalisation
##   - no clipping
##   - min contig 500bp
##   - RF mode
##   - non genome mode

echo
date
echo ======================= De Novo Assembly of Transcripts =======================

TRINITY_DIR=$OUT_DIR/Trinity
mkdir $TRINITY_DIR

Trinity \
 --seqType fq \
 --max_memory "$MEM"G \
 --min_contig_length 500 \
 --full_cleanup \
 --no_normalize_reads \
 --output $TRINITY_DIR/"$SAMPLE"_Trinity \
 --left $IG_READS_DIR/"$SAMPLE"_read1.fastq.gz \
 --right $IG_READS_DIR/"$SAMPLE"_read2.fastq.gz \
 --CPU $CORES

################################### Generate IG transcripts ##################################
## BLAST assembeled transcripts to the IMGT reference FASTA file
## Collect unique transcript IDs associated with IMGT reference

echo
date
echo ========================== IG Transcript extraction ===========================

for gene in ${IG_CHAIN[@]}
do

echo Extracting $gene transcripts:

if [[ "$gene" == "IGH" ]]; then
    IMGT_REF="$IGH_REF"
elif [[ "$gene" == "IGKL" ]]; then
    IMGT_REF="$IGKL_REF"
fi

blastn -db $IMGT_REF -query $TRINITY_DIR/"$SAMPLE"_Trinity.Trinity.fasta -outfmt 6 | \
  cut -f1 | \
  uniq | \
  xargs -n 1 samtools faidx $TRINITY_DIR/"$SAMPLE"_Trinity.Trinity.fasta > $OUT_DIR/"$SAMPLE"_"$gene"_transcripts.fasta

count=$(grep ">" $OUT_DIR/"$SAMPLE"_"$gene"_transcripts.fasta | wc -l)

echo "  - $count Transcripts Extracted"
echo

done

###################################### Clean Trinity Files #####################################
if ! $KEEP_INT; then
    if [[ -e "$TRINITY_DIR" ]]; then
        rm -r $TRINITY_DIR
    fi
fi

echo
echo ==============================================================================
echo "IgSeqR - QUANTIFY, ANNOTATE AND FILTER TRANSCRIPTS"
echo ==============================================================================

################################# Quantify IGH transcripts #################################
## Quantify extracted IG associated transcripts using Kallisto pseudoalignment

echo
date
echo =========================== Kallisto Quantification ===========================

for gene in ${IG_CHAIN[@]}
do
echo "Running Kallisto for $gene Transcripts..."
KALLSTO_DIR=$OUT_DIR/kallisto/$gene
mkdir -p $KALLSTO_DIR

## Build kallisto index for transcripts
kallisto index -i $KALLSTO_DIR/"$gene"_transcripts.index $OUT_DIR/"$SAMPLE"_"$gene"_transcripts.fasta

## Quantify number of reads mapping to transcripts
kallisto quant \
  -i $KALLSTO_DIR/"$gene"_transcripts.index \
  -o $KALLSTO_DIR \
  -t $CORES \
  $IG_READS_DIR/"$SAMPLE"_read1.fastq.gz $IG_READS_DIR/"$SAMPLE"_read2.fastq.gz

done

###################################### Clean IG Read Files #####################################
if ! $KEEP_INT; then
    if [[ -e "$IG_READS_DIR" ]]; then
        rm -r $IG_READS_DIR
    fi
fi
  
#for gene in "${IG_CHAIN[@]}"
#do
#  KALLSTO_DIR=$OUT_DIR/kallisto/$gene
#  if [ -d "$KALLSTO_DIR" ]; then
#    tail -n +2 "$KALLSTO_DIR/abundance.tsv" | sort -t $'\t' -k5,5nr | head -5 | cut -f1 | xargs -I {} samtools faidx "$OUT_DIR/${SAMPLE}_${gene}_transcripts.fasta" {} > "$OUT_DIR/${SAMPLE}_${gene}_TPM_filtered.fasta"
#    mv "$KALLSTO_DIR/abundance.tsv" "$OUT_DIR/${SAMPLE}_${gene}_abundance.tsv"
#    echo
#  else
#    echo "ERROR: The directory '$KALLSTO_DIR' does not exist."
#    exit 1
#  fi
#done
echo
date
echo ========================== Processing IgSeqR Outputs ==========================

for gene in "${IG_CHAIN[@]}"
do
  KALLSTO_DIR=$OUT_DIR/kallisto/$gene
  if [ -d "$KALLSTO_DIR" ]; then
    echo "Processing gene: $gene ..."
    
    # Define the FASTA file
    fasta_file="$OUT_DIR/${SAMPLE}_${gene}_transcripts.fasta"
    if [ ! -f "$fasta_file" ]; then
      echo "ERROR: FASTA file '$fasta_file' not found."
      continue
    fi
    
    # Collect all gene transcripts and merge their abundances into a single file
    report_file="$OUT_DIR/${SAMPLE}_${gene}_report.tsv"
    echo -e "sample_id\ttranscript_id\tlength\teff_length\test_counts\tTPM\tsequence" > "$report_file"
    echo " - Generating IgSeqR report and FASTA files"
    
    # Extract sequences for all transcripts and merge with abundance data
    if [ -f "$KALLSTO_DIR/abundance.tsv" ]; then
      while IFS=$'\t' read -r transcript_id length eff_length est_counts TPM; do
        sequence=$(samtools faidx "$fasta_file" "$transcript_id" | grep -v "^>" | tr -d '\n')
        echo -e "$SAMPLE\t$transcript_id\t$length\t$eff_length\t$est_counts\t$TPM\t$sequence" >> "$report_file"
      done < <(tail -n +2 "$KALLSTO_DIR/abundance.tsv")
    else
      echo "ERROR: Abundance file '$KALLSTO_DIR/abundance.tsv' not found."
      continue
    fi
    
    # Create a new FASTA file using the report
    modified_fasta="$OUT_DIR/${SAMPLE}_${gene}_transcripts.fasta"
    > "$modified_fasta"
    while IFS=$'\t' read -r sample_id transcript_id length eff_length est_counts TPM sequence; do
      header=">${sample_id}|${transcript_id}"
      echo -e "$header\n$sequence" >> "$modified_fasta"
    done < <(tail -n +2 "$report_file")
    
    # Filter for the top 5 most abundant transcripts by TPM
    dominant_report_file="$OUT_DIR/${SAMPLE}_${gene}_dominant_report.tsv"
    echo -e "sample_id\ttranscript_id\tlength\teff_length\test_counts\tTPM\tsequence" > "$dominant_report_file"
    echo " - Generating dominant IgSeqR report and FASTA from top 5 transcripts by TPM"
    tail -n +2 "$report_file" | sort -t $'\t' -k6,6nr | head -5 >> "$dominant_report_file"
    
    # Create a new FASTA file using the dominant report
    tpm_filtered_fasta="$OUT_DIR/${SAMPLE}_${gene}_TPM_filtered.fasta"
    > "$tpm_filtered_fasta"
    while IFS=$'\t' read -r sample_id transcript_id length eff_length est_counts TPM sequence; do
      header=">${sample_id}|${transcript_id}"
      echo -e "$header\n$sequence" >> "$tpm_filtered_fasta"
    done < <(tail -n +2 "$dominant_report_file")
    
  else
    echo "ERROR: The directory '$KALLSTO_DIR' does not exist."
    exit 1
  fi
done

###################################### Clean Kallisto Files #####################################

if ! $KEEP_INT; then
    if [[ -e "$OUT_DIR/kallisto/" ]]; then
        rm -r $OUT_DIR/kallisto/
    fi
fi

echo 
echo ==============================================================================
echo IgSeqR COMPLETE
echo ==============================================================================
echo