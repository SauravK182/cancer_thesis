#!/bin/bash

#########################
# Author: Saurav Kiri
# Last modified: 2022-12-20
# Description: Calls ATAC-seq peaks with MACS2
# Dependencies: MACS2, MultiQC (for summarizing MACS2 peak calls)
# Type bash call_atac_peaks.sh --h (or --help) to see list of all options
#########################

arg0="bash $(basename "$0" .sh).sh"   # (basename arg suffix) removes leading directory of arg and trailing suffix
blnk=$(echo "$arg0" | sed 's/./ /g')

# $* means space separated string of all arguments
error() {
    echo
    echo "ERROR $arg0: $*" >&2
    echo
    exit 1
}

usage_info() {
    echo "Usage: $arg0 [--help] --input {FILE 1 FILE 2 ... | INPUT_DIR}"
    echo "       $blnk --outdir OUTFILE_DIR --reportdir REPORTS_DIR [--single]"
    echo "       $blnk [--files] [--qval VALUE] [--broadthreshold VALUE]"
}

usage() {
    usage_info
    echo
    echo "For detailed help, please type:"
    echo "  bash $arg0 --help"
    exit 1
}

help() {
    usage_info
    echo
    echo "This script will automate the peak calling of ATAC-seq reads with MACS2, using --broad."
    echo "-------------------------------------------------------------------------------------------------"
    echo "There are other methods of calling ATAC-seq peaks, including using the options --nomodel"
    echo "--shift -100 --extsize 200. The shift -100 will shift the 5' ends first 100 bp in the 5' direction"
    echo "and extsize 200 will extend all reads to a fixed fragment size of 200 bp (5' to 3'). However, shift"
    echo "can only be set when BAMPE is OFF, and thus one can only do this analysis while in single-end mode."
    echo "The upside to this analysis is that due to the recentering of peaks, assuming the fragments are all 200"
    echo "bp (which is somewhat in line with that reported by Buenrostro et al., 2013), shifting 5' by 100 bp"
    echo "has the effect of re-centering fragments along the cutting sites. This script assumes the user is"
    echo "not interested in the cutting sites, but moreso in the paired-end fragment pileup."
    echo
    echo "Required arguments:"
    echo "  -i | --input           Pass input directory. With flag --files on, will assume a list of files for peak calling."
    echo "  -o | --outdir          Output directory for files produced by MACS2."
    echo "  -r | --reportdir       Output directory for MultiQC summary report."
    echo
    echo "Optional arguments:"
    echo "  -s | --single          Flag to indicate if ATAC-seq experiments were paired-end or not. Default = off (i.e., PE)."
    echo "  -f | --files           Flag indicating that --input parameter represents a list of files, not a directory. Default = off."
    echo "  -q | --qval            Optional argument to specify a q-value threshold for MACS2. Default is 0.05."
    echo "  -t | --broadthreshold  Optional argument for q-value threshold for broad peak calling. Default is 0.10."
    echo "  -h | --help            Print this message and exit."
    exit 0
}

# Set defaults
INPUTFILES=
OUTDIR=
REPORTDIR=
PAIRED='true'
FILES='false'
QVAL=0.05
BROADVAL=0.10

flags() {
    while test $# -gt 0
    do
        case "$1" in
        (-i|--input)
            shift
            [ $# = 0 ] && error "No input files/directory specified"
            INPUTFILES=()
            while [[ "$1" != -* ]] && [ ! -z "$1" ]
            do
                INPUTFILES+=("$1")
                shift
            done;;
        (-o|--outdir)
            shift
            [ $# = 0 ] && error "No output directory specified"
            OUTDIR="$1"
            shift;;
        (-r|--reportdir)
            shift
            [ $# = 0 ] && error "No report directory specified"
            REPORTDIR="$1"
            shift;;
        (-s|--single)
            PAIRED='false'
            shift;;
        (-f|--files)
            FILES='true'
            shift;;
        (-q|--qval)
            shift
            [ $# = 0 ] && error "Must specify q-value when using -q."
            QVAL="$1"
            shift;;
        (-t|--broadthreshold)
            shift
            [ $# = 0 ] && error "Must specify broad threshold cutoff when using -t."
            BROADVAL="$1";;
        (-h|--help)
            help;;
        (*)
            usage;;
        esac
    done
}

# Get flags
flags "$@"

if [ -z $INPUTFILES ] || [ -z $OUTDIR ] || [ -z $REPORTDIR ]    # if any required args are missing, quit
then
    usage
fi

# If files flag is set to true, save file list
# Else, get the .bam files from the directory
if $FILES
then
    BAMFILES=${INPUTFILES[@]}
else
    BAMFILES=${INPUTFILES}*.bam
fi

# For each ATAC-seq file, call MACS2
# First, if single flag was given, set format to BAM; else, keep format BAMPE
$PAIRED && FORMAT=BAMPE         # if paired, set format = BAMPE
$(! $PAIRED) && FORMAT=BAM      # if not paired, set format = BAM (user specified --single)
for FILE in $BAMFILES
do
    PREFIX=$(basename $FILE .bam)
    macs2 callpeak -t $FILE -f $FORMAT -g hs -n $PREFIX --outdir $OUTDIR --keep-dup all --broad -q $QVAL --broad-cutoff $BROADVAL --bdg --SPMR  # note ATAC-seq has no controls
    sort -k1,1 -k2,2n ${OUTDIR}${PREFIX}_treat_pileup.bdg > ${OUTDIR}${PREFIX}_treat_sorted.bdg
    bedGraphToBigWig ${OUTDIR}${PREFIX}_treat_sorted.bdg ~/project/hg19-build/contig_sizes.txt ${OUTDIR}${PREFIX}.bw
done

# Create MultiQC report of number of peaks called by MACS2 for each run
multiqc --module macs2 --force --outdir $REPORTDIR -n "multiqc_macs2_atac_peaks" $OUTDIR