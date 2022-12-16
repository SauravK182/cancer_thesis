#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-12
# Description: Align downloaded single- or paired-end RNA-seq FASTQ files to hg19 genome using HISAT2
# Dependencies: HISAT2, samtools (to produce sorted bam), RNASeQC (for quality checking), MultiQC (for report aggregation)
# Arguments: 
    # 1) path to directory with all FASTQ files 
    # 2) Path/basename of HISAT2 index 
    # 3) Path to .bam output directory
    # 4) Path to collapsed gtf (for RNASeQC)
    # 5) path to report outputdir 
    # 6) flag (se or pe) to indicate whether files are single-end or paired-end. Assumes paired-end by default.
# Usage: bash ~/project/scripts/rna-seq/rnaseq_align.sh <path-to-fastq> <fa-reference> <path-to-bam-output> <path-to-collapsed-gtf> <path-to-report-dir> <pe/se>
#########################

# Get all paths from arguments
READS_DIR=${1}
REF_PATH=${2}
BAM_DIR=${3}
ANNO_DIR=${4}
QC_DIR=${5}
FLAG="${6:-pe}"


# Set directory to that with FASTQ reads
cd ${READS_DIR}
echo "The current directory of FASTQ files is: ${PWD}"

# If flag is specified as SE, get FASTQ filenames and align with the file
# Else, get prefixes and align with both forward and reverse files

if [ ${FLAG} == "se" ]
then
    PREFIX_LIST=$(ls *fastq* | sed -r "s/[.]fastq.*//")
    for prefix in ${PREFIX_LIST}
    do
        FULLNAME=$(ls ${prefix}*)
        echo -e "The current alignment is taking place on: ${FULLNAME}\n"
        
        # Use -U to indicate unpaired read
        hisat2 --new-summary --summary-file "${BAM_DIR}${prefix}.txt" -p 6 \
        -x ${REF_PATH} -U "${FULLNAME}" | samtools sort -o "${BAM_DIR}${prefix}.bam"
    done
else
    # First, obtain the list of all paired-end read files
    ## These files should be of the form <base-name>_1.fastq and <base-name>_2.fastq per experiment
    ## sed is a text stream editor - can search and edit files without opening the file in a text editor
    ## sed usage is as follows for replacing strings: sed s/<instance-to-rep>/<replacement-str>; -r indicates use extended regexp
    ## uniq filters/removes ADJACENT duplicate lines from a text file
    ## Command due to https://www.biostars.org/p/98222/
    PREFIX_LIST=$(ls *fastq* | sed -r "s/_[12][.]fastq.*|_[12]_trimmed[.]fastq.*//" | uniq)

    # Therefore above, we are getting the list of all fastq files, replacing the "_1.fastq" or "_2.fastq" with nothing, then removing all duplicate names
    # This gives us only the prefixes we want - e.g., hpne-r1, etc.

    echo -e "The sequencing experiments to be aligned are: ${PREFIX_LIST}\n"

    # Align to genome with HISAT2
    for prefix in ${PREFIX_LIST}
    do
        READ_FILES=($(ls ${prefix}*))
        echo -e "The current alignment is taking place on: ${READ_FILES[0]} and ${READ_FILES[1]}\n"

        hisat2 --new-summary --summary-file "${BAM_DIR}${prefix}.txt" -p 6 \
        -x ${REF_PATH} -1 ${READ_FILES[0]} -2 ${READ_FILES[1]} | samtools sort -o "${BAM_DIR}${prefix}.bam"
    done
fi


# Perform RNASeQC quality check
for prefix in ${PREFIX_LIST}
do
    echo -e "Now running RNASeQC quality check on: ${prefix}\n"
    rnaseqc --sample ${prefix} ${ANNO_DIR} "${BAM_DIR}${prefix}.bam" ${QC_DIR}
done

# Aggregate reports from RNASeQC and summary HISAT2 files with MultiQC
echo -e "Completed alignment and RNASeQC. Aggregating reports with MultiQC\n"
multiqc --module hisat2 --module rna_seqc -o ${QC_DIR} -n "multiqc_alignment_results" \
-c ~/multiqc_configs/multiqc_config.yaml ${BAM_DIR} ${QC_DIR}
