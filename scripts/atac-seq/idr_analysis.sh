#!/bin/bash

#########################
# Author: Saurav Kiri
# Date last modified: 2023-04-09
# Description: Performs IDR analysis to keep only reproducible peaks between experiments
# Dependencies: IDR (https://github.com/nboley/idr), bedtools
# Type bash idr_analysis.sh --h (or --help) to see list of all options
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
    echo "Usage: $arg0 [-f] [-i IDR_VAL] -o -t <TREATMENT_1, TREATMENT_2...> -c CONTROL"
    echo "       $blnk [-a FASTA_FILE] POS_ARG"
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
    echo "Performs IDR analysis on sorted narrowPeak files from MACS2 to merge replicates and keep only"
    echo "peaks with a high-confidence of being reproducible (i.e., true binding sites). Bedtools is also"
    echo "used to find peaks that are specific to a treatment (or any condition) relative to the control"
    echo "condition. Finally, a FASTA file can be specified with -a to map identified reproducible peaks"
    echo "to sequences with bedtools"
    echo "***********************************************************************************************"
    echo "It should be noted that the simple intersection of peaksets to determine unique peaks is"
    echo "not generally recommended, but in this case, I wished to look for peaks that are specific to"
    echo "one treatment over a control - representing de novo accessible regions."
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo "  -f|--files          Indicates positional argument POS_ARG is a list of files."
    echo "  -i|--idr_val        IDR cutoff for determining reproducible peaks. Default is 0.05."
    echo "  -o|--outdir         Output directory for IDR results."
    echo "  -t|--treatment      Space-separated list of names representing the treatments. Note that EACH"
    echo "                      associated .narrowPeak treatment files should have this string in its name"
    echo "                      to be matched exactly. Program will FAIL if this does not occur."
    echo
    echo "  -c|--control        String representing the control files. Note that the associated .narrowPeak"
    echo "                      files should have this string in its filename AS AN EXACT MATCH. Currently only"
    echo "                      supports a single control at the moment."
    echo
    echo "  -a|--fasta          Optional genome FASTA file to be provided to map reproducible peaks to sequences."
    echo "  POS_ARG             Either a directory of .narrowPeak files (assumes they are sorted and end in _sorted.narrowPeak)"
    echo "                      or a list of _sorted.narrowPeak files (if flag -f is specified)."
    exit 0
}

containsElement () {
  local e match="$1"        # set local vars e (empty) and match
  shift                     # shift so only array elements remain in positional parameters
  for e; do [[ "$e" == "$match" ]] && return 0; done    # omitting 'in' within for loop will cause loop to loop over positional params
  return 1
}

# Set variable defaults/initialize variables
ISFILES=false
IDR_VAL=0.05
OUTDIR=
TREATMENTS=()
CONTROL=
FASTA=

# Create flags function to grab passed parameters
# See here for skeleton: https://stackoverflow.com/questions/192249/
flags() {
    # Create POS_ARGS variable to take in positional arguments
    POS_ARGS=()
    while test $# -gt 0
    do
        case "$1" in
        (-f|--files)
            ISFILES=true
            shift;;
        (-i|--idr_val)
            shift
            IDR_VAL="$1"
            shift;;
        (-o|--outdir)
            shift
            OUTDIR="$1"
            shift;;
        (-t|--treatment)
            shift
            while [[ "$1" != -* ]] && [[ ! -z "$1" ]]; do
                TREATMENTS+=("$1")
                shift
            done;;
        (-c|--control)
            shift
            CONTROL="$1"
            shift;;
        (-a|--fasta)
            shift
            FASTA="$1"
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

# Quit if outdir, treatments, or control is not specified
if [[ -z $OUTDIR ]] || [[ -z $TREATMENTS ]] || [[ -z $CONTROL ]]; then
    error "A necessary parameter was not specified. Please consult the documentation (flag -h) for details on usage."
fi

# Get treatment files
if ! $ISFILES; then
    PEAKFILES=${POS_ARGS}*_sorted.narrowPeak
else
    PEAKFILES=${POS_ARGS[@]}
fi

# Perform reverse sort on treatments
# This ensures that treatments with the same base name (e.g., POU and POU51) are not double counted later
ALL_NAMES=$( echo "${TREATMENTS[@]}" "$CONTROL" | sort -r )

# Create hash table to store replicates with associated string
declare -A PAIRS
for NAME in $ALL_NAMES; do
    FILES_MATCHING=()
    for FILE in $PEAKFILES; do
    # check for exact match b/w filename and names passed and see if file is already present in the array
        if [[ "$FILE" =~ $NAME ]] && ! containsElement $FILE ${PAIRS[@]}; then
            FILES_MATCHING+=($FILE)
        fi
    done
    [ ! -z $FILES_MATCHING ] && PAIRS+=([$NAME]=${FILES_MATCHING[@]})
done

# Run IDR to get reproducible peaks - capture stderr output as log
for KEY in ${!PAIRS[@]}; do
    echo "Now running IDR on samples for $KEY."
    idr --samples ${PAIRS[$KEY]} \
    --input-file-type narrowPeak \
    --rank q.value \
    --output-file "${OUTDIR}${KEY}" \
    -i $IDR_VAL \
    --plot 2> "${OUTDIR}${KEY}.log"
    # Extract first 3 BED cols and peak summit (10th col)
    cut -f 1,2,3,10 "${OUTDIR}${KEY}" > "${OUTDIR}${KEY}.bed"
done

# Run bedtools intersect to get peaks specific to the treatments, grab 200 bp window around summit for motif analysis
for TREAT in "${TREATMENTS[@]}"; do
    bedtools intersect -v -a "${OUTDIR}${TREAT}.bed" -b "${OUTDIR}${CONTROL}.bed" | awk '{ print $1, $2+$4-100, $2+$4+100 }' OFS='\t' > "${OUTDIR}${TREAT}_unique.bed"
done

# If user supplied genome FASTA file, get the associated DNA sequences for the unique peaks
if [ -z $FASTA ]; then
    echo "FASTA file not provided. Skipping peak FASTA generation."
elif [ -f $FASTA ]; then
    echo "Generating FASTA files..."
    for NAME in ${TREATMENTS[@]}; do
        bedtools getfasta -fi $FASTA -bed "${OUTDIR}${NAME}_unique.bed" -fo "${OUTDIR}${NAME}_sequences.fasta"
    done
else
    echo "Cannot access passed FASTA file."
fi
