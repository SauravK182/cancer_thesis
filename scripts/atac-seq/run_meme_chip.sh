#!/bin/bash

#########################
# Author: Saurav Kiri
# Date last modified: 2023-04
# Description: Runs MEME-ChIP and TomTom on identified enrichment sites to discover TF binding motifs.
# Dependencies: MEME Suite Software, TeXLive epstopdf
# Type bash run_meme_chip.sh --h (or --help) to see list of all options
#########################


arg0="bash $(basename "$0")"
blnk=$(echo "$arg0" | sed 's/./ /g')

error() {
    echo
    echo "ERROR $arg0: $*" >&2
    echo
    exit 1
}

usage_info() {
    echo "Usage: $arg0 -f FASTA -d DATABASE -o OUTPUTDIR [-w]"
}

usage() {
    echo
    echo "For detailed help, please type:"
    echo "	bash $arg0 --help"
    exit 1
}

help() {
    usage_info
    echo
    echo "Simple script for automating a run of MEME-ChIP followed by TomTom to discover enriched"
    echo "biological TF motifs in a given FASTA file. Note this script assumes the following:"
    echo "1) FASTA file is over the DNA alphabet"
    echo "2) You wish to call E-values (proportion of motifs with equal or greater log-odds ratio given a background) at 0.05"
    echo "3) You want to use the di-nucleotide permuted background rather than a normal background"
    echo "4) You have a single database file (e.g., a JASPAR .meme file)."
    echo "Script also produces .eps logo plots and converts the .eps to .pdf files with TeXLive's epstopdf."
    echo 
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo "-f | --fasta          FASTA file for testing motif enrichment."
    echo "-d | --database       Database file for finding matches to TF binding sites."
    echo "-o | --outputdir      Output directory. TomTom results will be in a subdirectory 'TomTom' here."
    echo "-w | --overwrite      Whether or not to overwrite the output directory if it exists. Default = off."
    echo "                      Note if -w is OFF, then the directory specified in -o should not exist."
    exit 0
}

# Set variable defaults/initialize variables
FASTA=
DATABSE=
OUTPUTDIR=
OVERWRITE=false

# Create flags function to grab passed parameters
# See here for skeleton: https://stackoverflow.com/questions/192249/
flags() {
    while test $# -gt 0
    do
        case "$1" in
        (-f|--fasta)
            shift
            FASTA="$1"
            shift;;
        (-d|--database)
            shift
            DATABASE="$1"
            shift;;
        (-o|--outputdir)
            shift
            OUTPUTDIR="$1"
            shift;;
        (-w|--overwrite)
            shift
            OVERWRITE=true
            shift;;
        (-h|--help)
            help;;
        (-*|--*)
            usage;;
        (*)
            shift;;
        esac
    done
}

# Get input variables using the flags function
flags "$@"

# Quit if one of the variables was not defined
[ -z $FASTA ] || [ -z $DATABASE ] || [ -z $OUTPUTDIR ] && error "One of the necessary variables was not set. Please see documentation (-h) for assistance."

# Run MEME-ChIP and TomTom
if $OVERWRITE; then DIRTYPE=-oc; else DIRTYPE=-o; fi
echo "Starting MEME-ChIP on $FASTA.."
meme-chip $DIRTYPE $OUTPUTDIR $FASTA && \
echo "Done. Starting TomTom analysis." && \
tomtom -o ${OUTPUTDIR}/tomtom/ ${OUTPUTDIR}/combined.meme $DATABASE -eps && \
for FILE in ${OUTPUTDIR}/tomtom/*.eps; do epstopdf $FILE; done      # Convert eps to pdf