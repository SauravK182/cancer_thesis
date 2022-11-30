#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-08
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

# Use FastQC on each FASTQ file in the directory
## Use REPORTS_DIR to direct the QC reports to user-defined dir
for f in ${FASTQ_DIR}*.fastq
do
    echo "Processing ${f}"
    fastqc --outdir ${REPORTS_DIR} ${f}
    echo -e "Finished QC on ${f}.\n The report is placed in ${REPORTS_DIR}."
done


# Use MultiQC to aggregate FastQC reports into a single report
## --force will force overwrite files
## Specify --module fastqc to tell multiqc this is only FastQC data
multiqc --module fastqc --force --outdir ${REPORTS_DIR} -n "multiqc_fastqc_results" ${REPORTS_DIR}
