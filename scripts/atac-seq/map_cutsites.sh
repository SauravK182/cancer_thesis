#!/bin/bash

#########################
# Author: Saurav Kiri
# Date last modified: 2023-04-09
# Description: Maps Tn5 transposase cut sites on MACS2 using narrow peak calling
# Dependencies: MACS2
# Type bash map_cutsites.sh --h (or --help) to see list of all options
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
    echo "Usage: $arg0 {-d || -f} -o OUTDIR [-q QVAL] [-s SHIFT] [-e EXTSIZE] POS_ARG"
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
    echo "Script designed for mapping Tn5 cutsites with ATAC-seq by ignoring model building, smoothening signal, and calling narrow peaks."
    echo "One may also be interested in using the -f BAMPE mode and allowing MACS to compute peaks with its default behavior. However, I wrote"
    echo "this script envisioning mapping cut sites for motif analysis to determine what, if any, regulatory elements are becoming de novo accessible"
    echo "across a treatment."
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo "  -d | --dir          Indicates positional argument POS_ARG is a directory of .bam files."
    echo "  -f | --files        Indicates positional argument POS_ARG is a list of .bam files."
    echo "  -o | --outdir       Specifies output directory for peak files."
    echo "  -q | --qval         FDR-adjusted cutoff for peak significance. Default is 0.10."
    echo "  -s | --shiftsize    Parameter to pass to MACS2 --shift for moving back reads. Default is -100 bp."
    echo "  -e | --extsize      Parameter for smoothing read density signal with MACS2 param --extsize. Default is 200 bp."
    exit 0
}

# Set variable defaults/initialize variables
ISDIR=false
ISFILES=false
OUTDIR=
QVAL=0.1
SHIFT=-100
EXTSIZE=200


# Create flags function to grab passed parameters
# See here for skeleton: https://stackoverflow.com/questions/192249/
flags() {
    # Create POS_ARGS variable to take in positional arguments
    POS_ARGS=()
    while test $# -gt 0; do
        case "$1" in
        (-d|--dir)
            ISDIR=true
            shift;;
        (-f|--files)
            ISFILES=true
            shift;;
        (-o|--outdir)
            shift
            OUTDIR="$1"
            shift;;
        (-q|--qval)
            shift
            QVAL="$1"
            shift;;
        (-s|--shift)
            shift
            SHIFT="$1"
            shift;;
        (-e|--extsize)
            shift
            EXTSIZE="$1"
            shift;;
        (-h|--help)
            help;;
        (-*|--*)
            usage;;
        (*)
            # In this case, we have reached positional args; save them
            POS_ARGS+=("$1")
            shift;;
        esac
    done
}

# Get input variables using the flags function
flags "$@"

# Stop if user did not specify output dir and used neither or both of -d and -f
if [[ $ISDIR == $ISFILES ]] || [[ -z $OUTDIR ]]; then
    error "Improper argument list. Please see the documentation via -h or --help for details."
fi

# Get list of ATAC-seq bam files
if $ISFILES; then
    BAMFILES=${POS_ARGS[@]}
else
    BAMFILES=${POS_ARGS}*.bam
fi

# Call ATAC-seq peaks and sort narrowPeak file by -log10(q), which is 9th column
for FILE in $BAMFILES; do
    PREFIX=$(basename $FILE .bam)
    macs2 callpeak -t $FILE -f BAM -g hs -n $PREFIX --outdir $OUTDIR --keep-dup all --nomodel --shift $SHIFT --extsize $EXTSIZE -q $QVAL
    # sort file numerically by 9th col in descending order (largest to smallest)
    sort -k9,9nr "${OUTDIR}${PREFIX}_peaks.narrowPeak" > "${OUTDIR}${PREFIX}_peaks_sorted.narrowPeak"
done