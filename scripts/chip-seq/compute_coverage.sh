#!/bin/bash

#########################
# Author: Saurav Kiri
# Date last modified: 2023-04
# Description: Computes density matrices for ChIP/ATAC-seq data for downstream visualization.
# Dependencies: deepTools
# Type bash compute_coverage.sh --h (or --help) to see list of all options
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
    echo "Usage: $arg0 [-h] { -b BAMDIR | -w BIGWIGDIR } -o OUTDIR [-r BEDFILE] [-f OUTNAME]"
    echo "       $blnk [-a] [-n NORMALIZEUSING] [-s BINSIZE] [-l SMOOTHING] [-p]"
    echo "       $blnk [-t THREADS] [-c WINDOW] [-e EXTSIZE]"
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
    echo "Script designed to run bamCompare/bamCoverage & computeMatrix from deepTools in an automated fashion."
    echo "The computation of bamCompare or bamCoverage can be skipped if a directory of .bw files is provided."
    echo "Note this script assumes ChIP-seq data has input controls and ATAC-seq data does not. However, the -a"
    echo "flag can be used to specify whether or not there are input controls independent of experiment type."
    echo "In this script, it is assumed that the reference sequences for plotting (if desired) are around TSSs."
    echo 
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo "  -h | --help             Print this message and exit."
    echo "  -b | --bamdir           Directory of bam files. If -a not specified, assumes both exp and input controls are here"
    echo "  -w | --bigwigdir        Directory of .bw files for computeMatrix. Mutually exclusive with -b."
    echo "  -o | --outdir           Directory for output files of bamCoverage and/or computeMatrix."
    echo "  -r | --bedfile          BED file required for computeMatrix - if not passed, this will not be run."
    echo "  -f | --outname          Outfile name for computeMatrix file. Not required if -r not specified."
    echo "  -a | --atac             Indicates no input controls => bamCoverage is run instead of bamCompare. Default = off."
    echo "  -n | --normalization    Normalization strategy to use. Default is BPM."
    echo "  -s | --binsize          Bin size to use for read density calculation. Default is 50."
    echo "  -l | --smoothing        Window size to smooth read averages over. Default is 150."
    echo "  -p | --paired           Indicates .bam files are paired-end. Default is off."
    echo "  -t | --threads          Number of threads to run deepTools on. Default is 10."
    echo "  -c | --window           Window size for computing matrix. Default is 2000."
    echo "  -e | --extsize          Size for extending reads. Value ignored if -p is indicated. Default for single-end is 200."
    exit 0
}

containsElement () {
  local e match="$1"        # set local vars e (empty) and match
  shift                     # shift so only array elements remain in positional parameters
  for e; do [[ "$e" == "$match" ]] && return 0; done    # omitting 'in' within for loop will cause loop to loop over positional params
  return 1
}

# Function to separate a list of .bam files into a hashtable of input controls and associated treatments
parseBam() {
    local -a input treat
    local file test
    # since file is empty, for file loops through positional args
    for file; do
        if [[ $file = *input* ]]; then
            input+=($file)
        else
            treat+=($file)
        fi
    done

    declare -A -g PAIRS
    for file in ${input[@]}; do
        local trt=()
        local expr=$(basename $file .bam)
        local exprnum=${expr#*input}      # recall ${word#pattern} returns expression with shortest matching pattern deleted
        local exprtype=$(echo ${expr} | sed -r "s/-input.*|_input.*//")
        for test in ${treat[@]}; do
        # if the current treatment file is already a value in PAIRS, then skip
            if ! containsElement $test ${PAIRS[@]} && [[ "$test" = *$exprtype*$exprnum* ]]; then
                trt+=($test)
            fi
        done
        [ ! -z $trt ] && PAIRS+=([$file]=${trt[@]})  # append key-value (control-treatments) pair to associative array 
    done                                             # but only if inputs have corresponding treatments
}

# Set variable defaults/initialize variables
BAMDIR=
BWDIR=
OUTPUTDIR=
BEDFILE=
OUTNAME=
ATAC=false
NORMALIZE="BPM"
BINSIZE=50
SMOOTH=150
PAIRED=false
THREADS=10
WINDOW=2000
EXTSIZE=200

# Create flags function to grab passed parameters
# See here for skeleton: https://stackoverflow.com/questions/192249/
flags() {
    while test $# -gt 0; do
        case "$1" in
        (-b|--bamdir)
            shift
            BAMDIR="$1"
            shift;;
        (-w|--bigwigdir)
            shift
            BWDIR="$1"
            shift;;
        (-o|--outdir)
            shift
            OUTPUTDIR="$1"
            shift;;
        (-r|--bed)
            shift
            BEDFILE="$1"
            shift;;
        (-f|--outname)
            shift
            OUTNAME="$1"
            shift;;
        (-a|--atac)
            ATAC=true
            shift;;
        (-n|--normalization)
            shift
            NORMALIZE="$1"
            shift;;
        (-s|--binsize)
            shift
            BINSIZE="$1"
            shift;;
        (-l|--smoothing)
            shift
            SMOOTH="$1"
            shift;;
        (-p|--paired)
            PAIRED=true
            shift;;
        (-t|--threads)
            shift
            THREADS="$1"
            shift;;
        (-c|--window)
            shift
            WINDOW="$1"
            shift;;
        (-e|--extsize)
            shift
            EXTSIZE="$1"
            shift;;
        (-h|--help)
            help;;
        (*)
            usage;;
        esac
    done
}

# Get input variables using the flags function
flags "$@"

# Check to make sure bamdir and bwdir are not both set and all other variables are ready
if [ ! -z $BAMDIR ] && [ ! -z $BWDIR ]; then error "Must select one of bamdir or bwdir, but not both."; fi
[ -z $OUTPUTDIR ] && [ ! -z $BAMDIR ] && OUTPUTDIR=$BAMDIR
[ -z $OUTPUTDIR ] && [ ! -z $BWDIR ] && OUTPUTDIR=$BWDIR
if [ ! -z $BEDFILE ] && [ -z $OUTNAME ]; then error "Must specify an outfile name if a BED file is provided."; fi

# Set extendRead size - if paired, deepTools will automatically calculate extension size
if $PAIRED; then
    EXTSIZE=
elif ! $PAIRED && [ -z $EXTSIZE ]; then
    error "You did not specified .bam files are paired. In such a case, please indicate an extension size with -e."
fi


# If bam files are provided, get files and run bamCompare or bamCoverage
if [ ! -z $BAMDIR ]; then
    BAMFILES=$BAMDIR*.bam
    # if not atac (i.e., has inputs), get inputs and run bamCompare
    if ! $ATAC; then
        parseBam $BAMFILES
        # Run bamCompare on each treatment-input pair
        for CONTROL in "${!PAIRS[@]}"; do
            for TREAT in ${PAIRS[$CONTROL]}; do
            echo "Now running bamCompare on:"
            echo "$TREAT and $CONTROL"
            FILENAME=$( basename $TREAT .bam )
            # Run bamCompare with desired settings
            bamCompare -b1 "$TREAT" \
            -b2 "$CONTROL" \
            -o "${OUTPUTDIR}${FILENAME}.bw" \
            --binSize $BINSIZE \
            --normalizeUsing "$NORMALIZE" \
            --smoothLength $SMOOTH \
            --extendReads $EXTSIZE \
            --centerReads \
            -p $THREADS \
            --scaleFactorsMethod None 2> "${OUTPUTDIR}${FILENAME}.log" && echo -e "Done!\n"
            done
        done
    else
        for BAM in $BAMFILES; do
            FILENAME=$( basename $BAM .bam )
            echo "Now running bamCoverage on:"
            echo "$BAM"
            # Run bamCoverage for ATAC-seq data - note we do not extend reads since we only want to map the ends of the reads, guaranteed accessible
            bamCoverage -b $BAM \
            -o "${OUTPUTDIR}${FILENAME}.bw" \
            --binSize $BINSIZE \
            --normalizeUsing "$NORMALIZE" \
            -p $THREADS 2> "${OUTPUTDIR}${FILENAME}.log" && echo -e "Done!\n"
        done
    fi
fi

# If bedfile is passed, run computeMatrix
if [ ! -z $BEDFILE ]; then
    if [ -z $BWDIR ] && [ ! -z $BAMDIR ]; then
        BWDIR=$OUTPUTDIR
    fi
    BWFILES=$BWDIR*.bw
    echo "Now running computeMatrix with TSSs defined in $BEDFILE."
    computeMatrix reference-point --referencePoint TSS \
    -S $BWFILES \
    -b $WINDOW -a $WINDOW \
    -R $BEDFILE \
    --skipZeros \
    -o "${BWDIR}${OUTNAME}.gz" \
    -p $THREADS && echo -e "Done!\n"
fi