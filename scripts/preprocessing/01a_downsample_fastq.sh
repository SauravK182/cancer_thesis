#!/bin/bash

#########################
# Author: Saurav Kiri
# Date last modified: 2023-01-23
# Description: Downsamples FASTQ files to specified length
# Dependencies: Seqtk, pigz (if gzip desired)
# Type bash 01a_downsample_fastq.sh -h (or --help) to see list of all options
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
    echo "Usage: $arg0 [--help] [--files] [--single] [--raw]"
    echo "       $blnk --size [--outdir] {/path/to/files | file1 file2 ...}"
}

usage() {
    echo
    echo "For detailed help, please type:"
    echo "	$arg0 --help"
    exit 1
}

help() {
    usage_info
    echo
    echo "Uses the seqtk package to downsample a set of FASTQ files. Assumes users wishes to downsample to the same"
    echo "depth for all files, or at least to the same percentage of original depth for all files. Requires seqtk"
    echo "to be present in \$PATH. If --single is not specified, seeds will be automatically set such that paired"
    echo "reads are extracted from paired-end sequencing data."
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo "  -f | --files            Indicates positional argument(s) is/are a list of files instead of a directory."
    echo "  -s | --single           Indicates FASTQ files are single-end rather than paired-end."
    echo "                          Note that paired-end filenames are assumed to be of the form _1.fq, _2.fq, _1.fq.gz, or _2.fq.gz"
    echo "  -n | --size             Relative size or fraction of original reads to downsample to."
    echo "  -r | --raw              Indicates not to gzip output files. By default, pigz is used for gzipping."
    echo "  -o | --outdir           Output directory for downsampled FASTQ files."
    echo "                          If positional argument is a directory, this is optional, and by default, files will be written"
    echo "                          to the passed input directory. Required parameter if used with -f"
    echo "  -h | --help             Print this message and quit."
    echo "  positional              By default, a directory to search for FASTQ files. With -f, indicates a list of files."
    exit 0
}

# Set variable defaults/initialize variables
ISFILES='false'
SAMPLESIZE=
PAIRED='true'
GZIP='true'
OUTDIR=

# Create flags function to grab passed parameters
# See here for skeleton: https://stackoverflow.com/questions/192249/
flags() {
    # Create POS_ARGS variable to take in positional arguments
    POS_ARGS=()
    while test $# -gt 0
    do
        case "$1" in
        (-f|--files)
            ISFILES='true'
            shift;;
        (-n|--size)
            shift
            SAMPLESIZE="$1"
            shift;;
        (-s|--single)
            PAIRED='false'
            shift;;
        (-r|--raw)
            GZIP='false'
            shift;;
        (-o|--outdir)
            shift
            OUTDIR="$1"
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

# Get list of FASTQ file to be downsampled
if $ISFILES; then
    FASTQ_FILES=${POS_ARGS[@]}
else
    FASTQ_FILES=$(ls $POS_ARGS*fastq*)
    [ -z OUTDIR ] && OUTDIR=$POS_ARGS   # set default outdir if outdir was not specified by user
fi

echo -e "List of FASTQ files to be downsampled is:\n$FASTQ_FILES"

# Get number of unique seeds needed
if $PAIRED; then
    # Use wc -w to get number of words (space-separated elements) in FASTQ, divide by 2 bc PE
    NUM_SEEDS=$(($(wc -w <<< $FASTQ_FILES) / 2))
else
    NUM_SEEDS=$(wc -w <<< $FASTQ_FILES)
fi

# Create array of seeds, use shuf to create the random list of NUM_SEEDS elements
SEEDS=( $(shuf -i 0-10000 -n $NUM_SEEDS | sort -n) )

# Perform downsampling
[ ! -d "$OUTDIR/subsampled" ] && mkdir "$OUTDIR/subsampled" # make a sub-dir for subsampled files if it doesn't exist
echo -e "Writing sub-sampled files to $OUTDIR/subsampled.\n"
PREFIX_LIST=( $(echo -e "$FASTQ_FILES" | sed -r "s/_[12][.]fastq.*|[.]fastq.*//" | uniq) )
for INDEX in ${!PREFIX_LIST[@]}; do # ${!array[@]} returns list of valid indices
    if $PAIRED; then
        # Use printf '%s\n' to format entries; '%s' is placeholder for string input
        PAIREDFILES=( $(printf '%s\n' $FASTQ_FILES | grep ${PREFIX_LIST[$INDEX]}) ) # Search for matching files from prefix
        echo "Now running seqtk sample with seed ${SEEDS[$INDEX]} on:"
        echo "${PAIREDFILES[0]}"
        echo -e "${PAIREDFILES[1]}\n"
        # Use % parameter expansion to remove shortest matching pattern to rename files to "subsampled"
        BASEONE=$(basename ${PAIREDFILES[0]})
        BASETWO=$(basename ${PAIREDFILES[1]})
        seqtk sample -s ${SEEDS[$INDEX]} ${PAIREDFILES[0]} $SAMPLESIZE > "$OUTDIR/subsampled/${BASEONE%*_1.fastq*}_subsampled_1.fastq"
        seqtk sample -s ${SEEDS[$INDEX]} ${PAIREDFILES[1]} $SAMPLESIZE > "$OUTDIR/subsampled/${BASETWO%*_2.fastq*}_subsampled_2.fastq"
    else
        SINGLEFILE=$(printf '%s\n' $FASTQ_FILES | grep ${PREFIX_LIST[$INDEX]})
        echo "Now running seqtk sample with seed ${SEEDS[$INDEX]} on:"
        echo -e "$SINGLEFILE\n"
        seqtk sample -s ${SEEDS[$INDEX]} $SINGLEFILE $SAMPLESIZE > "$OUTDIR/subsampled/${SINGLEFILE%*.fastq*}_subsampled.fastq"
    fi
done

# # Run pigz if gzipping desired
if $GZIP; then
    for OUTFILE in $OUTDIR/subsampled/*subsampled*fastq; do
        pigz $OUTFILE
    done
fi