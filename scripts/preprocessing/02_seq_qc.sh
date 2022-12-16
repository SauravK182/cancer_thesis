#!/bin/bash

#########################
# Author: Saurav Kiri
# Last modified: 2022-12-14
# Description: Run FastQC on a directory of FASTQ files and use MultiQC to aggregate reports
# Dependencies: FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and MultiQC (https://multiqc.info)
# Arguments: 
    # 1) Directory with FASTQ files to perform quality check on
    # 2) directory to place QC reports in
# Usage: bash ~/project/scripts/preprocessing/02_rna_seq_fastqc.sh /path/to/fastqdir/ /path/to/reports/
#########################

# Get directory of FASTQ files and reports
FASTQ_DIR=${1}
REPORTS_DIR=${2}
echo "The directory of FASTQ files is set to: ${FASTQ_DIR}"

# Get list of trimmed FASTQ files and perform FastQC
# Set to run on 6 files simultaneously (threads = 6)
# Use two wildcards to match either .fastq.gz or .fastq
FASTQ_FILES=$(ls ${FASTQ_DIR}*.fastq*)
fastqc --threads 6 --outdir ${REPORTS_DIR} ${FASTQ_FILES}


# Use MultiQC to aggregate FastQC reports into a single report
## --force will force overwrite files
## Specify --module fastqc to tell multiqc this is only FastQC data
multiqc --module fastqc --force --outdir ${REPORTS_DIR} -n "multiqc_fastqc_results" ${REPORTS_DIR}
