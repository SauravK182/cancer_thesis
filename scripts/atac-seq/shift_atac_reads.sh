#!/bin/bash

#####################
# Author: Saurav Kiri
# Last modified: 2022-12-19
# Description: Uses alignmentSieve to perform +4/-5 bp shift for ATAC-seq reads
# Dependencies: deepTools
# Type bash shift_atac_reads.sh --help for more info
#####################

arg0="bash $(basename "$0" .sh).sh"   # (basename arg suffix) removes leading directory of arg and trailing suffix
blnk=$(echo "$arg0" | sed 's/./ /g')

error() {
    echo
    echo "ERROR $arg0: $*" >&2
    echo
    exit 1
}

usage_info() {
    echo "Usage: $arg0 {--indir INPUT_DIR | --infiles FILE 1 [FILE 2 ... FILE N]}"
    echo "       $blnk --outdir OUTPUT_DIR [--numthreads INT]"
}

usage() {
    usage_info
    echo
    echo "For detailed help, please type:"
    echo "  bash shift_atac_reads.sh --help"
    exit 1
}

help() {
    usage_info
    echo
    echo "Script for performing the 'conventional' +4 (forward)/-5 (reverse) bp ATAC-seq read shifting"
    echo "using --ATACshift alignmentSieve modeule from deepTools."
    echo "----------------------------------------------------------------------------------------------------"
    echo "It should be noted that in general, this +4/-5 bp shifting is done by convention from the original"
    echo "ATAC-seq paper (Buenrostro et al., 2013). The Tn5 transposase binds to DNA as a homodimer and makes"
    echo "nicks in the backbone approximately 9 bp apart, creating 9 bp ssDNA overhangs ('sticky ends') in the"
    echo "cut DNA, where it will then ligate the double-stranded mosaic-end adapters. Since the reads that will be"
    echo "sequencd from these adapter inserts will be from the actual 9 bp-offset cutting sites themselves, the"
    echo "authors of the original paper decided to shift reads mapping to the + strand by 4 bp downstream and"
    echo "reads mapping to the - strand by 5 bp upstream. The idea is to have the 5' ends of reads pileup at"
    echo "the center of a Tn5 cutting event. In truth, such a shift is likely to have a minimal effect on"
    echo "downstream analysis, but I decided to write this script to keep in line with convention."
    echo
    echo "Arguments (optional arguments in [brackets]):"
    echo "  -d | --indir         Input directory for .bam files to be extracted from."
    echo "  -f | --infiles       Alternatively, a set of .bam files can be specified directly."
    echo "  -o | --outdir        Where to place output .bam files."
    echo "  -n | --numthreads    Number of threads to run on. Default is 8."
    exit 0
}

INFILES=
INDIR=
OUTDIR=
NUMTHREADS=8

flags() {
    while test $# -gt 0
    do
        case "$1" in
        (-d|--indir)
            shift
            INDIR="$1"
            shift;;
        (-f|--infiles)
            shift
            INFILES=()
            while [[ "$1" != -* ]] && [ ! -z "$1" ] # loop over all listed files, breaking if we reach the end of passed parameters
            do
                INFILES+=("$1")
                shift
            done;;
        (-o|--outdir)
            shift
            OUTDIR="$1"
            shift;;
        (-n|--numthreads)
            shift
            NUMTHREADS="$1"
            shift;;
        (-h|--help)
            help;;
        (*)
            usage;;
        esac
    done
}

flags "$@"

# If passed args is files, save .bam files. Else, get .bam files from directory
if [ ! -z $INFILES ] && [ -z $INDIR ]
then
    BAMFILES=${INFILES[@]}
elif [ -z $INFILES ] && [ ! -z $INDIR ]
then
    BAMFILES=$(ls $INDIR*.bam)
else    # this implies the user set both infiles and indir or neither. Print error message and quit.
    error "Please specify either an input file(s) or a directory of .bam files (but not both)."
fi

# Run alignmentSieve to do ATAC-seq shift on each bam, save to outdir
# Get prefix, save for output file
for FILE in $BAMFILES
do
    PREFIX=$(basename $FILE .bam)   # basename gets rid of leading directory and suffix, will leave only file name
    if [ -f $FILE.bai ]             # alignmentSieve requires indexed .bam files - check if index exists before invoking samtools index
    then
        echo "File $FILE is already indexed. Moving on..."
    else
        echo -e "Now indexing $FILE for alignmentSieve.\n"
        samtools index -@ 6 $FILE
    fi
    echo -e "Now running alignmentSieve on $FILE with $NUMTHREADS threads and writing to ${OUTDIR}${PREFIX}-shifted.bam"
    alignmentSieve -b $FILE -o "${OUTDIR}${PREFIX}-shifted.bam" -p $NUMTHREADS --ATACshift
done
