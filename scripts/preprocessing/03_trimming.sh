#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-23
# Description: Trims reads in a passed directory for adapter sequences and quality
# Requires: BBDuk package from the BBTools suite
# Arguments: 
    # 1) directory to list of fastq files to trim
    # 2) output directory for trimmed files 
    # 3) q-score to trim to
    # 4) min length to filter by 
    # 5) directory for QC reports
    # 6) character string "se" or "pe" to indicate if files are single-end or paired end. Assumes paired-end by default.
# Usage: bash ~/project/scripts/preprocessing/03_trimming.sh <fastq-dir> <out-dir> <q-score> <minlen> <reports-dir> <se/pe>
#########################

# Get paths/scores from arguments
FASTQ_DIR=${1}
TRIMOUT_DIR=${2}
QSCORE=${3}
MINLEN=${4}
QC_DIR=${5}
FLAG="${6:-pe}" # Sets the 6th variable to "pe" by default

cd ${FASTQ_DIR}
echo -e "Current directory of FASTQ files is: ${FASTQ_DIR}"

# To account for possible .fastq or .fastq.gz extension, dynamically match extension
    # Get full file name
    # Use parameter expansion ${file#*.} to delete everything up to and including first period
    # Save this + leading period as file extension
# file_fullname will do this programmatically and will only create global variable EXT
# See https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html for more
file_fullname () {
    local FILE_PREFIX=${1}
    local FULLNAMES=($(ls ${FILE_PREFIX}*))
    EXT=()
    local name
    for name in ${FULLNAMES[@]}
    do
        EXT+=(".${name#*.}")
    done
}

# If files are SE (i.e., flag var is set to "se"), get list of SE FASTQ files
# Else, assumes reads are PE and trim/filter accordingly.
if [ ${FLAG} == "se" ]
then
    PREFIX_LIST=$(ls *fastq* | sed -r "s/[.]fastq.*//")
    echo -e "The FASTQ files to be trimmed are:\n${PREFIX_LIST}"

    # Use bbduk with vars to trim and filter reads; save stderr summary to file
    for prefix in ${PREFIX_LIST}
    do
        file_fullname ${prefix}
        # Returns an array containing the extensions for each file
        echo -e "Currently starting trim on FASTQ file for ${prefix}${EXT[0]}.\n"

        bbduk.sh in="${prefix}${EXT[0]}" out="${TRIMOUT_DIR}${prefix}_trimmed${EXT[0]}" \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        ref=/home/sk_ounce/bbmap/resources/adapters.fa \
        qtrim=rl trimq=${QSCORE} \
        minlen=${MINLEN} 2> "${QC_DIR}${prefix}_trimstats.txt"

        echo -e "\nDone with trim on FASTQ file for ${prefix}. The trimmed file is in ${TRIMOUT_DIR}.\n"
    done
else
    PREFIX_LIST=$(ls *fastq* | sed -r "s/_[12][.]fastq.*//" | uniq)
    echo -e "The FASTQ files to be trimmed are:\n${PREFIX_LIST}"

    # Use bbduk with vars to trim and filter reads; save stderr summary to file
    for prefix in ${PREFIX_LIST}
    do
        file_fullname ${prefix}
        # In this case, the array returned should have 2 elements
        # Since for each PE experiment (prefix), there is a forward and reverse read file
        echo -e "Currently starting trim on FASTQ files for ${prefix}.\n"

        bbduk.sh in1="${prefix}_1${EXT[0]}" in2="${prefix}_2${EXT[1]}" \
        out1="${TRIMOUT_DIR}${prefix}_1_trimmed${EXT[0]}" out2="${TRIMOUT_DIR}${prefix}_2_trimmed${EXT[1]}" \
        ktrim=r \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe tbo \
        ref=/home/sk_ounce/bbmap/resources/adapters.fa \
        qtrim=rl \
        trimq=${QSCORE} \
        minlen=${MINLEN} 2> "${QC_DIR}${prefix}_trimstats.txt"

        echo -e "\nDone with trim on FASTQ files for ${prefix}. The trimmed files are in ${TRIMOUT_DIR}.\n"
    done
fi


# Get list of trimmed FASTQ files and perform FastQC
# Set to run on 6 files simultaneously (threads = 6)
# Use two wildcards to match either .fastq.gz or .fastq
fastqc --threads 6 --outdir ${QC_DIR} ${TRIMOUT_DIR}*.fastq*

# Generate aggregate MultiQC report
multiqc --module fastqc --module bbmap --force --outdir ${QC_DIR} \
-c ~/multiqc_configs/multiqc_config.yaml -n "multiqc_trimmed_results" ${QC_DIR}