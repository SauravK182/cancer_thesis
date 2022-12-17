#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-12-15
# Description: Calls ChIP-seq peaks with MACS2
# Dependencies: MACS2, MultiQC (for summarizing MACS2 peak calls)
# Arguments: 
    # 1) -f OR -d to indicate whether you are passing a directory or set of .bam files
    # 2) Set of .bam files (-f) OR directory for .bam files (-d)
    # 3) bins for multiBamSummary
    # 4) Directory to place coverage matrix from multiBamSummary and for plotCorrelations output
    # 5) Name of matrix output + heatmap pdf
# Usage: bash ~/project/scripts/rna-seq/rnaseq_align.sh [-f <bam-files> | -d <bam-dir>] <bins> <deeptools-outdir> <matrix-pdf-name>
#########################