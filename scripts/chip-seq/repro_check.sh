#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-12-14
# Description: Checks reproducibility of read coverage using plotCorrelation from deepTools
# Dependencies: samtools (for indexing files), deepTools (for using plotCorrelations)
# Arguments: 
    # 1) -f OR -d to indicate whether you are passing a directory or set of .bam files
    # 2) Set of .bam files (-f) OR directory for .bam files (-d)
    # 3) bins for multiBamSummary
    # 4) Directory to place coverage matrix from multiBamSummary and for plotCorrelations output
    # 5) Name of matrix output + heatmap pdf
# Usage: bash ~/project/scripts/rna-seq/rnaseq_align.sh [-f <bam-files> | -d <bam-dir>] <bins> <deeptools-outdir> <matrix-pdf-name>
#########################

# Initialize variables to avoid contamination from env vars
IS_FILES='false'
IS_DIR='false'

# Get the flag arguments if set
while getopts 'fd' flag; do # Note: use a colon (f:) if a flag is expecting an argument
    case "${flag}" in       # case is used to simply nested if
    f) IS_FILES='true' ;;   # of the form pattern) statements ;;
    d) IS_DIR='true' ;;     # if pattern matches, set IS_DIR to 'true'
    esac
done

# Assign BAM_FILES based on passed flag
if ${IS_FILES}
then
    BAM_FILES=${2}
elif ${IS_DIR}
then
    BAM_DIR=${2}
    BAM_FILES=${BAM_DIR}*.bam
else    # Quit if user does not specify
    echo "Please specify whether you are passing a directory of .bam files (-d) or a list of bam files (-f)."
    exit 1
fi

# Get rest of variables from input
BINS=${3}
OUTDIR=${4}
NAME=${5}

# Index all bam files
for file in ${BAM_FILES}
do
    if [ -f ${file}.bai ] # if index already exists (-f checks if file exists in regular format), skip
    then
        echo -e "${file} is already indexed. Moving on to the next file...\n"
    else                # else, index the file
        echo -e "Now indexing ${file}.\n"
        samtools index -@ 6 ${file}
    fi
done

# Call multiBamSummary
echo -e "Now running multiBamSummary\n"
multiBamSummary bins --binSize ${3} --bamfiles ${BAM_FILES} -o "${OUTDIR}${NAME}.npz"

# Use produced .npz to make a correlation plot
echo -e "Finished multiBamSummary. Now creating heatmap from ${OUTDIR}${NAME}.npz"
plotCorrelation --corData "${OUTDIR}${NAME}.npz" \
--corMethod pearson \
--whatToPlot heatmap \
--plotFile "${OUTDIR}${NAME}.pdf" \
--skipZeros \
