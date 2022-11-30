#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-13
# Description: Downloads and runs QC on ChIP-seq data corresponding to input, CTCF, and K27ac from Ren et al. (2021) pancreatic cancer SRA, GSE149103
# Requires: Script files fastq_download.sh (downloads FASTQ) and seq_dq.sh (runs FastQC/MultiQC quality check)
# Usage: bash ~/project/scripts/chip-seq/ren_download_qc.sh
#########################

# Set path to FASTQ dir - used as args by both scripts
CHIP_DIR=/mnt/d/SK/data/ren-panc/chip-seq/

# Specify file within this dir for SRR names + file names for download
CHIP_NAMES=chip-k27ac-ctcf-input-acc.csv

# Set paths for FastQC/MultiQC output
OUT_DIR=~/project/QC/Ren_2021/ChIP-seq/Pre-align/

# Call both scripts to download and run QC
bash ~/project/scripts/preprocessing/01_fastq_download.sh ${CHIP_DIR} ${CHIP_NAMES}
bash ~/project/scripts/preprocessing/02_seq_qc.sh ${CHIP_DIR} ${OUT_DIR}