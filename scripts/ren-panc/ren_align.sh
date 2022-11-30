#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-12
# Description: Uses written scripts rnaseq_align.sh and bam_to_counts.sh to align and perform QC of Ren et al RNA-seq data
# Usage: bash ~/project/scripts/ren-panc/ren_align.sh
#########################

# Assign variables to pass to rnaseq_align.sh
FASTQ=~/project/data/ren-panc/rna-seq/raw/
REFERENCE=~/project/hg19-build/grch37_tran/genome_tran
OUT_BAM=/mnt/d/SK/data/ren-panc/rna-seq/aligned/
ANNOTATION=~/project/hg19-build/human_hg19_annotations_collapsed.gtf
QC=~/project/QC/Ren_2021/RNA-seq/Post-align/

# Assign variables to pass to bam_to_counts.sh
## Note that OUT_BAM above acts as our BAM folder and output dir
FC_FILENAME="ren_rna_rawcounts.txt"
GTF_INPUT=~/project/hg19-build/Homo_sapiens.GRCh37.87.gtf
FC_QC_DIR=~/project/QC/Ren_2021/RNA-seq/Post-featurecounts/

# Call scripts with assigned variables
bash ~/project/scripts/rna-seq/rnaseq_align.sh ${FASTQ} ${REFERENCE} ${OUT_BAM} ${ANNOTATION} ${QC}
bash ~/project/scripts/rna-seq/bam_to_counts.sh ${OUT_BAM} ${OUT_BAM} ${FC_FILENAME} ${GTF_INPUT} ${FC_QC_DIR}

# Aggregate all QC reports into one