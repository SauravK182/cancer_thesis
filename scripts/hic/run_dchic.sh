#!/bin/bash

#########################
# Author: Saurav Kiri
# Date last modified: 2/24/23
# Description: Will run dcHiC on a selected directory of HiC-Pro matrix files
# Dependencies: dcHiC, remove_Y.r and blacklist_Y.r Rscripts
# Type bash run_dchic.sh --h (or --help) to see list of all options
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
    echo "Usage: $arg0 -d HICPRO_DIR -r RES -o DCHIC_OUTDIR -s DCHIC_SCRIPT_DIR"
    echo "       $blnk [-g GENOME] [-i INPUTNAME] [-b] [-y] [-h]"
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
    echo "Script for automating the preparation of .matrix and .bed files from HiC-Pro for dcHiC as well as calling dcHiC itself."
    echo "This script assumes the HiC folder names are of the form <sample-name>-hic-<replicate-number>. Additionally, it is"
    echo "assumed that the HiC-Pro directory hiearachy is relatively unchanged and that the .matrix file of interest is named exactly as"
    echo "HiC-Pro named it, i.e., <full-sample-name>_<res>.matrix. The same assumption is made for the .bed files."
    echo "This script automatically produces the input file required by dcHiC using the sample names; this .txt file is"
    echo "tab-separated and will be placed in the output directory specified by the -o parameter."
    echo "Currently, this script will only call the --pcatype options of cis, select, analyze, subcomp, and viz (in that order) for dcHiC."
    echo
    echo "Arguments (optional arguments in [brackets] above):"
    echo -e "  -d, --hicdir        Top-level directory for HiC-Pro; this folder should contain the hic_results folder from HiC-Pro.\n"
    echo -e "  -r, --res           Resolution to call differential compartments at; will be used to choose the appropriate .matrix/.bed files.\n"
    echo -e "  -o, --hicout        Output directory for where to place the dcHiC differential call results.\n"
    echo -e "  -s, --script        Directory in which the dchicf.r script is located."
    echo -e "  -g, --genome        Genome build for the given HiC-Pro style files. Used by dcHiC to correlate PC1/PC2 to GC content and gene density"
    echo -e "                      to select the PC which properly represents A/B compartments. Default = hg19.\n"
    echo -e "  -i, --inputfile     Name of input file to be created for dcHiC. Will automatically overwrite file with the same name."
    echo -e "                      Default = 'dchic_input.txt'.\n"
    echo -e "  -b, --blacklist     Flag for whether or not to run the blacklist_Y.r script to blacklist the Y chromosome, as well as"
    echo -e "                      add the prefix 'chr' to all chromosome names in the .bed file, required for dcHiC. Default = FALSE.\n"
    echo -e "  -y, --rmY           Flag for whether or not to run the remove_Y.r script to remove the Y chromosome from matrix and bed files."
    echo -e "                      May help for preventing errors during --pcatype cis in dcHiC. Default = FALSE."
    exit 0
}

# Set variable defaults/initialize variables
TOPDIR=
RES=
DCOUTDIR=
DCSCRIPT=
GENOME="hg19"
INPUTNAME="dchic_input.txt"
BLACKLIST='false'
RMY='false'

# Create flags function to grab passed parameters
# See here for skeleton: https://stackoverflow.com/questions/192249/
flags() {
    while test $# -gt 0; do
        case "$1" in
        (-d|--hicdir)
            shift
            TOPDIR="$1"
            shift;;
        (-r|--res)
            shift
            RES="$1"
            shift;;
        (-o|--hicout)
            shift
            DCOUTDIR="$1"
            shift;;
        (-s|--script)
            shift
            DCSCRIPT=$1
            shift;;
        (-g|--genome)
            shift
            GENOME="$1"
            shift;;
        (-i|--inputfile)
            shift
            INPUTNAME="$1"
            shift;;
        (-b|--blacklist)
            BLACKLIST='true'
            shift;;
        (-y|--rmY)
            RMY='true'
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

# If any required flags are missing, quit
if [ -z $TOPDIR ] || [ -z $RES ] || [ -z $DCOUTDIR ] || [ -z $DCSCRIPT ]; then
    error "Missing one or more parameters. Type $arg0 -h to see the documentation."
fi

# Run blacklisting/removing Y if necessary, create input sheet for dcHiC
> "$DCOUTDIR/$INPUTNAME"
MATRIXDIR="$TOPDIR/hic_results/matrix"
SAMPLES=( $(ls $MATRIXDIR) )
for EXP in ${SAMPLES[@]}; do
    DATADIR="$MATRIXDIR/$EXP/raw/$RES/"
    MATRIXFILE="$DATADIR/${EXP}_$RES.matrix"
    BEDFILE="$DATADIR/${EXP}_${RES}_abs.bed"

    echo "Current matrix and bed files are:"
    echo "$MATRIXFILE"
    echo -e "$BEDFILE\n"

    if $BLACKLIST; then
        echo -e "Now running blacklist_Y.r on $BEDFILE.\n"
        Rscript ~/project/scripts/hic/blacklist_Y.r $BEDFILE
        BEDFILE=$(echo $BEDFILE | sed "s/.bed/_blacklisted.bed/")
        echo $BEDFILE
    fi

    if $RMY; then
        echo "Now running remove_Y.r on:"
        echo "$MATRIXFILE"
        echo -e "$BEDFILE\n"
        Rscript ~/project/scripts/hic/remove_Y.r -m $MATRIXFILE -b $BEDFILE
        MATRIXFILE=$(echo $MATRIXFILE | sed -r "s/.matrix$/_filtered.matrix/")
        BEDFILE=$(echo $BEDFILE | sed "s/.bed/_filtered.bed/")
    fi

    # Get replicate number, experiment name; do some defense to make sure exp name is valid
    EXPNAME=$( echo $EXP | sed -r "s/-hic-[0-9]+$//" | tr '-' '_' | tr '.' '_' )
    if [[ $EXPNAME =~ ^[0-9] ]]; then EXPNAME="Sample_$EXPNAME"; fi # append sample_ if exp name starts w/ digit
    REPNUM=${EXP##*-}   # assuming user defined samples properly
    REPNAME="${EXPNAME}_R$REPNUM"

    # Append to input file and tell user
    echo -e "$MATRIXFILE\t$BEDFILE\t$REPNAME\t$EXPNAME" >> "$DCOUTDIR/$INPUTNAME"
    echo "The metadata for $EXP has been added to $INPUTNAME, including:"
    echo "Matrix file: $MATRIXFILE"
    echo "BED file: $BEDFILE"
    echo "Replicate name: $REPNAME"
    echo -e "Sample name: $EXPNAME\n"
done


# Run dcHiC
cd "$DCOUTDIR"
echo -e "Now running dcHiC...\n"

Rscript $DCSCRIPT/dchicf.r --file $INPUTNAME --pcatype cis; echo -e "Done with --pcatype cis. Moving on to --pcatype select...\n"
Rscript $DCSCRIPT/dchicf.r --file $INPUTNAME --pcatype select --genome $GENOME; echo -e "Done with --pcatype select. Moving on to --pcatype analyze...\n"
Rscript $DCSCRIPT/dchicf.r --file $INPUTNAME --pcatype analyze; echo -e "Done with --pcatype analyze. Moving on to --pcatype subcomp...\n"
Rscript $DCSCRIPT/dchicf.r --file $INPUTNAME --pcatype subcomp; echo -e "Done with --pcatype subcomp. Moving on to --pcatype viz...\n"
Rscript $DCSCRIPT/dchicf.r --file $INPUTNAME --pcatype viz --genome $GENOME
echo -e "Done!\n"