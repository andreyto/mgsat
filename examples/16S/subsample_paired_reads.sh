#!/bin/bash
#' Subsample all FASTQ files with the same random seed
#' This ensures that paired end reads stay paired
#' @param seed Random seed to use
#' @param read_count count or ratio of reads to take from each file
#' as used by `seqtk sample`
#'
#' Input file names are read from stdin
#' Output files are created in cwd with the same basename as inputs
#' This script uses BASH features
seed=$1
shift
[ -n "$seed" ] || exit 1
read_count=$1
shift
[ -n "$read_count" ] || exit 1

while read file_inp; do
    file_out=$(basename "$file_inp")
    if ! [ -e "$file_inp" ]; then
        (>&2 echo "Input file '$file_inp' does not exist")
    fi
    if [ "$file_inp" -ef "$file_out" ]; then
        (>&2 echo "Output file name $file_out points to input file name $file_inp")
    fi
    echo "Subsampling $file_inp"
    seqtk sample -s"$seed" "$file_inp" "$read_count" > "$file_out"
done

