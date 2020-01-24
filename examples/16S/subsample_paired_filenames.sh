#!/bin/bash
#' Subsample paired file names (FASTQ files) by taking matched pairs
#' @param seed Random seed to use
#' @param count max count of file pairs to output
#'
#' Input file names are read from stdin
#' Output file name are written to stdout
seed=$1
shift
[ -n "$seed" ] || exit 1
count=$1
shift
[ -n "$count" ] || exit 1

sort | xargs -n2 | shuf | head -n "$count" | tr '[:blank:]' '\n'

