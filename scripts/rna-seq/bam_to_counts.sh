#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-12
# Description: Uses featureCounts to create final raw count matrix from aligned .bam data (assumes paired-end)
# Dependencies: featureCounts and MultiQC; requires also a valid list of .bam files
# Arguments:
    # 1) directory of bam files
    # 2) directory for output of featureCounts
    # 3) filename for featurecounts
    # 4) path to .gtf annotation file
    # 5) path to desired multiqc output dir
    # 6) String "se" or "pe", indicating whether reads are single-end or paired-end. Assumes paired-end by default.
# Usage: bash ~/project/scripts/rna-seq/bam_to_counts.sh <bam-dir> <output-dir> "output-filename.txt" <path-to-gtf> <multiqc-output-dir> <se/pe>
#########################

# Get dirs/filename from arguments
INPUT_DIR=${1}
OUT_DIR=${2}
OUTNAME=${3}
GTF_FILE=${4}
QC_OUT=${5}
FLAG="${6:-pe}"

# Get all .bam files from directory
cd ${INPUT_DIR}
echo -e "Searching for the list of .bam files in: ${PWD}\n"
BAM_LIST=$(ls *.bam)

# Call featureCounts with list of .bam files
    # If files are paired-end, specify -p to count fragments instead of reads
    # Else, call featureCounts without -p to indicate single-end reads
if [ ${FLAG} == "pe" ]
then
    featureCounts -p -a ${GTF_FILE} -o "${OUT_DIR}${OUTNAME}" ${BAM_LIST}
else
    featureCounts -a ${GTF_FILE} -o "${OUT_DIR}${OUTNAME}" ${BAM_LIST}
fi

# Run MultiQC to aggregate report from .summary files
multiqc --module featureCounts -o ${QC_OUT} -n "multiqc_featurecounts_results" ${OUT_DIR}
