#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-08
# Description: Downloads and runs QC on all RNA-seq data from Ren et al. (2021) pancreatic cancer SRA, GSE149103
# Requires: Script files fastq_download.sh (downloads FASTQ) and seq_dq.sh (runs FastQC/MultiQC quality check)
# Usage: bash ~/project/scripts/preprocessing/ren_download_qc.sh
#########################

# Set path to FASTQ dir, which is used as arg for both scripts
RNA_DIR=~/project/data/ren-panc/rna-seq/

# Specify file with SRR names for download
SRR_FILE=ren_rnaseq_acclist.txt

# Set path to reports dir, used as arg for seq_qc.sh
OUTPUT_DIR=/mnt/c/Users/jvons/Documents/NCF/Thesis/QC/Ren_2021/RNA-seq/

# Call both scripts to download, then run QC on FASTQ data
bash fastq_download.sh ${RNA_DIR} ${SRR_FILE}
bash seq_qc.sh  ${RNA_DIR} ${OUTPUT_DIR}
