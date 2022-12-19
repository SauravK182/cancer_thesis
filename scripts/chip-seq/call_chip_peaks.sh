#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-12-18
# Description: Calls ChIP-seq peaks with MACS2
# Dependencies: MACS2, MultiQC (for summarizing MACS2 peak calls)
# Type bash call_chip_peaks.sh --h (or --help) to see list of all options
#########################


arg0=$(basename "$0" .sh)   # Removes leading directory and trailing suffix
blnk=$(echo "$arg0" | sed 's/./ /g')

# $* means space separated string of all arguments
error() {
    echo
    echo "ERROR $arg0: $*" >&2
    echo
    exit 1
}

usage_info() {
    echo "Usage: $arg0 [--help] --input {FILE 1 FILE 2 ... | INPUT_DIR} [--control FILE 1 FILE 2 ...]"
    echo "       $blnk --outdir OUTFILE_DIR --reportdir REPORTS_DIR {--broad | --narrow} [--files]"
    echo "       $blnk [--qval VALUE] [--broadthreshold VALUE]"
}

usage() {
    usage_info
    echo
    echo "For detailed help, please type:"
    echo "  bash call_chip_peaks.sh --help"
    exit 1
}

help() {
    usage_info
    echo
    echo "Script for automating peak calling with MACS2. Assumes the use of model or -f BAMPE (with flag --paired)."
    echo "This script will attempt to automatically call peaks for treatment pairs and their corresponding input control."
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo "  -i | --input           Pass input directory. With flag --files on, will assume a list of files for peak calling."
    echo "                         Assumes the input control file(s) are included with the list/dir if --control not specified."
    echo
    echo "  -c | --control         Input control file(s). If not specified, will assume the input controls are in --input."
    echo "                         Note this script assumes all input control files have 'input' in their name."
    echo
    echo "  -o | --outdir          Output directory for files produced by MACS2."
    echo "  -r | --reportdir       Output directory for MultiQC summary reports."
    echo "  -b | --broad           Indicates to use option --broad when peak calling with MACS2. Mutually eclusive with --narrow."
    echo "  -n | --narrow          Indicates to use narrow peak calling with MACS2. Mutually exclusive with --broad."
    echo "  -f | --files           Flag indicating that --input parameter represents a list of files, not a directory. Default = off."
    echo "  -q | --qval            Optional argument to specify a q-value threshold for MACS2. Default is 0.05."
    echo "  -t | --broadthreshold  Optional argument for threshold for broad peak calling (only for --broad). Default is 0.10."
    echo "  -h | --help            Print this message and exit."
    exit 0
}

# Set variable defaults/initialize variables
INPUT=
CONTROL=
ISFILE='false'      # Assume passed input is directory by default
ISBROAD='false'
ISNARROW='false'    # set narrow and broad to false to initialize
QVAL=0.05
CUTOFF=0.1
PAIRED='false'      # assume single-end by default
OUTDIR=
REPORTDIR=

# Skeleton from: https://stackoverflow.com/questions/64257286/

# Best to wrap vars in quotes when they can be empty or have spaces
# Shell can break strings with spaces into multiple args if not quoted
# Use case to check for cases of the first positional argument and save flag value (if applicable)
# Is basically short-hand for a long if/else if chain
flags() {
    while test $# -gt 0     # while number of passed parameters is greater than 0
    do
        case "$1" in        # will repeatedly loop the first positional argument through these patterns to find a match - almost like else if
        (-i|--input)
            shift   # consume the flag
            [ $# = 0 ] && error "No input files/directory specified"    # && only does cmd 2 if cmd 1 was successful
            INPUT=()
            while [[ "$1" != -* ]]      # in order to capture multiple files, keep adding filenames until we reach next flag
            do                          # Next flag starts with -, use globstar -* to indicate we have reached the next parameter
                INPUT+=("$1")           # See http://mywiki.wooledge.org/BashFAQ/031 for [ ] vs. [[ ]]
                shift
            done;;          # ;; tells case statement to stop searching if it finds a match
        (-c|--control)
            shift
            [ $# = 0 ] && error "Must specify input control files if using -c"
            CONTROL=()
            while [[ "$1" != -* ]]
            do
                CONTROL+=("$1")
                shift
            done;;
        (-f|--files)
            ISFILE='true'
            shift;;
        (-b|--broad)
            ISBROAD='true'
            shift;;
        (-n|--narrow)
            ISNARROW='true'
            shift;;
        (-q|--qval)
            shift
            [ $# = 0 ] && error "Must specify a q-value if using -q"
            QVAL="$1"
            shift;;
        (-t|--broadthreshold)
            shift
            [ $# = 0 ] && error "Must specify a broad cutoff if using -t"
            CUTOFF="$1"
            shift;;
        (-p|--paired)
            PAIRED='true'
            shift;;
        (-o|--outdir)
            shift
            [ $# = 0 ] && error "Must specify an output directory for MACS2 peak calling files"
            OUTDIR="$1"
            shift;;
        (-r|--reportdir)
            shift
            [ $# = 0 ] && error "Must specify an output directory for MultiQC summary report"
            REPORTDIR="$1"
            shift;;
        (-h|--help)
            help;;
        (*)     # for anything else, print the usage statement and quit with exit status 1
            usage;;
        esac
    done
}

flags "$@"      # call flags with argument list

# Use || for "or" in bash - will only execute second command if first never completes
# -z checks if a variable is unset or set to empty string ""
# If any of the required variables are unset (user did not specify), print usage and quit
if [ -z $INPUT ] || [ -z $OUTDIR ] || [ -z $REPORTDIR ]
then
    usage
fi

# Check to make sure broad and narrow are not the same (i.e., user did not specify or did not set)
if [ $ISBROAD = $ISNARROW ]
then
    error "Please specify EITHER --broad or --narrow. The options are mutually exclusive."
fi

# List all parameters for user (may remove in future iterations)
echo
echo "input is ${INPUT[@]}"
echo "isfile is $ISFILE"
echo "isbroad is $ISBROAD"
echo "isnarrow is $ISNARROW"
echo "qval is $QVAL"
echo "cutoff is $CUTOFF"
echo "paired is $PAIRED"
echo "outdir is $OUTDIR"
echo "reportdir is $REPORTDIR"
echo

sleep 3s    # pause for 3s to let user review the parameters

# If input is a list of files, save them, else, get a list of all .bam files from the directory
# Note that if [ ], when given a string to test (which has alias [ ]), will check if string is non-empty
# If string is non-empty, it will do the command
# To use true/false literals, simply use if true or if false; avoid using test
if $ISFILE
then
    BAM_FILES=${INPUT[@]}
else
    BAM_FILES=${INPUT}*.bam
fi

# Get input control file(s) and separate from treatment files
CONTROL=()
TREATMENT=()
for FILE in $BAM_FILES
do
    if [[ $FILE == *input*]]    # if "input" is in the name of the file, add file to the input control array
    then
        CONTROL+=($FILE)
    else
        TREATMENT+=($FILE)
    fi
done

# For each file in BAM_FILES, call MACS2 and save output to outdir with prefix name
for FILE in $BAM_FILES
do
    PREFIX=$(echo $FILE | sed -r "s/[.]bam//")
    if $PAIRED
    then
        if $ISBROAD
        then
            echo "Now running MACS2 with options -f BAMPE, --broad, and -q $QVAL and --broad-cutoff $CUTOFF"
            macs2 callpeak 