#!/bin/bash

#########################
# Author: Saurav Kiri
# Date: 2022-11-10
# Description: Converts .txt files edited in Windows to a UNIX-friendly version
    ## Lines of windows .txt files end with "\r\n", where "\r" represents a carriage return
    ## However on UNIX, .txt files only end with the newline character "\n"
# Arguments: Input name of/path to text file to be converted; output name/path to output file
# Usage: bash ~/project/scripts/preprocessing/win_to_nix_file.sh /path/to/winfile /path/to/outputfile

# Grab input file and output name from args
WIN_FILE=${1}
NIX_FILE=${2}

# Use the text manipulation command awk
## "\r$" is regexp - note that $ will match the position before the FIRST newline in the string
## Therefore, we are effectively telling sub() to match the \r before the \n at the end of the string
## And replace it with nothing
## Then we print the output and write it all to a new file

awk '{ sub("\r$", ""); print }' ${WIN_FILE} > ${NIX_FILE}