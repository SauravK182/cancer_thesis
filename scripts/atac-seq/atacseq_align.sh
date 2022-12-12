#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-12-11
# Description: Align downloaded single- or paired-end ATAC-seq FASTQ files to hg19 genome using BWA-MEM
# Dependencies: BWA-MEM, samtools (to produce sorted bam), MultiQC (for report aggregation)
# Arguments: 
    # 1) path to directory with all FASTQ files 
    # 2) path/basename of bwa index files 
    # 3) path to .bam output directory
    # 4) path to directory for MultiQC report to be placed
    # 5) flag (se or pe) to indicate whether files are single-end or paired-end. Assumes paired-end by default.
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

# Get files corresponding to blacklisted regions
BLACKLIST=$(ls ~/project/blacklist/*reformatted.bed)

# Generate sequence of canonical nuclear chromosomes to keep in final .bam
CHROM=($(seq 1 1 22) X Y)

# Get prefix list (filenames without extensions)
# Since some files might be named, e.g., _1.fastq or _trimmed.fastq.gz, need to account for this
PREFIX_LIST=$(ls *fastq* | sed -r "s/[.]fastq.*|_[12][.]fastq.*|_[12]_trimmed[.]fastq.*//" | uniq)

# If flagged specified to SE, get FASTQ names and do unpaired alignment
# Else, assume PE (default), get prefixes and align with both forward and reverse reads
for prefix in ${PREFIX_LIST}
do
    if [ ${FLAG} == "se" ]
    then
        # Accounts for if file is .fastq or .fastq.gz
        FULLNAME=$(ls ${prefix}*)
        echo -e "The current alignment is taking place on: ${FULLNAME}\n"

        # Align, filter via pipes
            # First, remove secondary alns, unampped alns, and alns w/ MAPQ < 20
            # Then remove reads mapping to consensus blacklist regions
        # Use tee to create an intermediate file for getting base mapping stats
        bwa mem -t 6 ${GENOME_INDEX} ${FULLNAME} | samtools sort --threads 6 | \
        tee "${BAM_OUTDIR}${prefix}-intermediate.bam" | \
        samtools view -@ 6 -F 260 -q 20 -b - | \
        bedtools intersect -v -abam stdin -b ${BLACKLIST} > "${BAM_OUTDIR}${prefix}-prefil.bam"
    else
        READ_FILES=($(ls ${prefix}*))
        echo -e "The current alignment is taking place on: ${READ_FILES[0]} and ${READ_FILES[1]}\n"
        
        # Use forward and reverse reads to map to genome
        # Align, filter via pipes
            # First, remove secondary alns, unampped alns, and alns w/ MAPQ < 20
            # Keep only reads that are properly paired
            # Then remove reads mapping to consensus blacklist regions
        # Use tee to create an intermediate file for getting base mapping stats
        bwa mem -t 6 ${GENOME_INDEX} ${READ_FILES[0]} ${READ_FILES[1]} | samtools sort --threads 6 | \
        tee "${BAM_OUTDIR}${prefix}-intermediate.bam" | \
        samtools view -@ 6 -F 260 -f 0x02 -q 20 -b - | \
        bedtools intersect -v -abam stdin -b ${BLACKLIST} > "${BAM_OUTDIR}${prefix}-prefil.bam"
    fi
    # Use samtools to keep only nuclear chromosomes
    # Note that this step requires an indexed bam file
    samtools index "${BAM_OUTDIR}${prefix}-prefil.bam"
    samtools view -@ 6 -b "${BAM_OUTDIR}${prefix}-prefil.bam" ${CHROM[@]} > "${BAM_OUTDIR}${prefix}.bam"

    # Get mapping statistics with samtools stats/idxstats and final reads kept w/ flagstat
    # Note MultiQC requires idxstats in filename to recognize idxstats
    samtools stats "${BAM_OUTDIR}${prefix}-intermediate.bam" > "${BAM_OUTDIR}${prefix}.txt"
    samtools idxstats "${BAM_OUTDIR}${prefix}-prefil.bam" > "${BAM_OUTDIR}${prefix}-idxstats.txt"
    samtools flagstat "${BAM_OUTDIR}${prefix}.bam" > "${BAM_OUTDIR}${prefix}-filtered.txt"

    rm "${BAM_OUTDIR}${prefix}-intermediate.bam"
    rm ${BAM_OUTDIR}${prefix}-prefil*     # Removes .bam and .bai
done

# Aggregate reports from samtools stats into one with MultiQC
multiqc --force --module samtools -o ${REPORTS_DIR} -n "multiqc_alignment_results" ${BAM_OUTDIR}