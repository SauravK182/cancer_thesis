#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-07
# Description: Download the SRA FASTQ files for a given set of SRRs and outfile names
# Dependencies: SRA Toolkit (https://github.com/ncbi/sra-tools/wiki)
# Arguments: 
    # 1) A .csv with SRR filename and output file names for fasterq-dump
        # Note: ensure this .csv does not possess the Windows-specific carriage return (\r) newline char 
    # 2) Directory to place dumped FASTQ-files
    # 3) (Optional) "gzip" - will indicate to gzip the fastq files
# Usage: bash ~/project/scripts/preprocessing/01_fastq_download.sh /path/to/download/dir/ "filename.csv"
#########################

# Get file and directory path from argument
DATA_DIR=${1}
FILE_NAME=${2}
GZIP_FLAG=${3:-"null"}

# Set directory to dir for sequence download from argument
cd ${DATA_DIR}
echo -e "The current data directory is: ${PWD}\n" # -e means "enables interpretation of backslash escapes"

# Read in file to get SRR and matching outfile name
## Note IFS is "internal field separator", determines how to do word splitting
## -r: do not allow backslashes to escape characters
## The final "<" redirects a file into the command
## E.g., sort < words.txt will apply sort to words.txt
while IFS=, read -r srr exp
do
    # prefetch.2.8.0 had certification failure, use prefetch.3.0.0 instead to pre-cache SRR files
    echo -e "Currently working on ${srr}, with outfile name prefix ${exp}\n"
    prefetch.3.0.0 --progress ${srr}

    # Download FASTQ from fetched file
    fasterq-dump.3.0.0 --progress --force --outfile "${exp}.fastq" ${srr}
    echo -e "Downloaded FASTQ files for ${srr} (outfile name: ${exp})\n"

    # Delete folder containing .sra file since it is no longer required
    rm -r ${srr}
done < ${FILE_NAME}

# Gzip files using parallel gzip (pigz) if requested
if [ ${GZIP_FLAG} == "gzip"]
then
    for file in ${DATA_DIR}*.fastq
    do
        pigz ${file}
    done
fi