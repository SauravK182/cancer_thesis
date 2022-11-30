#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-23
# Description: Align downloaded single- or paired-end ChIP-seq FASTQ files to hg19 genome using BWA-MEM
# Dependencies: BWA-MEM, samtools (to produce sorted bam), MultiQC (for report aggregation)
# Arguments: 
    # 1) path to directory with all FASTQ files 
    # 2) path/basename of bwa index files 
    # 3) path to .bam output directory 
    # 4) flag (se or pe) to indicate whether files are single-end or paired-end. Assumes paired-end by default.
# Usage: bash ~/project/scripts/rna-seq/rnaseq_align.sh <path-to-fastq> <name-of-fa-reference> <path-to-bam-output> <path-to-report-dir> <se/pe>
#########################

# Get filenames/paths from arguments
FASTQ_DIR=${1}
GENOME_INDEX=${2}
BAM_OUTDIR=${3}
REPORTS_DIR=${4}
FLAG="${5:-pe}"

# Set directory to reads w/ FASTQ files
cd ${FASTQ_DIR}
echo -e "The current directory of FASTQ files is: ${PWD}"

# If flagged specified to SE, get FASTQ names and do unpaired alignment
# Else, assume PE (default), get prefixes and align with both forward and reverse reads
if [ ${FLAG} == "se" ]
then
    PREFIX_LIST=$(ls *fastq | sed -r "s/[.]fastq//")
    for prefix in ${PREFIX_LIST}
    do
        echo -e "The current alignment is taking place on: ${prefix}.fastq\n"

        # Align, pipe to samtools to sort
        bwa mem -t 6 ${GENOME_INDEX} "${prefix}.fastq" | samtools sort --threads 6 -o "${BAM_OUTDIR}${prefix}-intermediate.bam"

        # Get mapping statistics with samtools stat
        samtools stats "${BAM_OUTDIR}${prefix}-intermediate.bam" > "${BAM_OUTDIR}${prefix}.txt"

        # Filter .bam file and remove the intermediate file
            # Remove secondary alignments, unmapped reads and alignments w/ MAPQ < 20
        samtools view -@ 6 -F 260 -q 20 -b "${BAM_OUTDIR}${prefix}-intermediate.bam" > "${BAM_OUTDIR}${prefix}.bam"
        samtools flagstat "${BAM_OUTDIR}${prefix}.bam" > "${BAM_OUTDIR}${prefix}.txt"
        rm "${BAM_OUTDIR}${prefix}-intermediate.bam"
    done
else
    # Since some files might be named, e.g., _1.fastq or _1_trimmed.fastq, need to account for this
    PREFIX_LIST=$(ls *fastq | sed -r "s/_[12][.]fastq|_[12]_trimmed[.]fastq//" | uniq)

    # Each filehandle should be unique to the forward and reverse reads
    # Obtain all files with handle = prefix and place into an array (indexing starts at 0)
    for prefix in ${PREFIX_LIST}
    do
        READ_FILES=($(ls ${prefix}*))
        echo -e "The current alignment is taking place on: ${READ_FILES[0]} and ${READ_FILES[1]}\n"
        
        # Use forward and reverse reads to map to genome
        bwa mem -t 6 ${GENOME_INDEX} ${READ_FILES[0]} ${READ_FILES[1]} | samtools sort --threads 6 -o "${BAM_OUTDIR}${prefix}-intermediate.bam"

        # Get mapping statistics w/ samtools stat
        samtools stats "${BAM_OUTDIR}${prefix}-intermediate.bam" > "${BAM_OUTDIR}${prefix}.txt"

        # Filter .bam file and remove the intermediate file
            # Remove secondary alignments, unmapped reads and alignments w/ MAPQ < 20
            # Keep only reads aligned as part of proper pair
        samtools view -@ 6 -F 260 -f 0x02 -q 20 -b "${BAM_OUTDIR}${prefix}-intermediate.bam" > "${BAM_OUTDIR}${prefix}.bam"
        samtools flagstat "${BAM_OUTDIR}${prefix}.bam" > "${BAM_OUTDIR}${prefix}-filtered.txt"
        rm "${BAM_OUTDIR}${prefix}-intermediate.bam"
    done
fi

# Aggregate reports from samtools stats into one with MultiQC
multiqc --module samtools -o ${REPORTS_DIR} -n "multiqc_alignment_results" ${BAM_OUTDIR}