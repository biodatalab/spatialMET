#!/bin/bash

##Use this (as example) to run in cluster
# Cluster pixels (i.e., spatial domain/niche detection) from intensity data using hcdist
# Input is a text file containing a table of peak intensities
# Second argument is a string specifying the name of output file containing 
# the pixel tree.
# Third argument specifies the number of cores to use.
#
# go_detect_spatial_domains.sh /home/user01/intensity_table.tsv_filtered.tsv /home/user01/intensity_table_feature_tree.tre 32
#

# Check if the file path argument is provided
#re='^[0-9]+$'
#if [ -z "$1" ] || [[ $1 =~ $re ]]; then
#	echo "Invalid input path. Usage: $0 <file_path> <out_path> <cores>"
#	exit 1
#fi

# Check if the output file path s provided
#if [ -z "$2" ] || [[ $2 =~ $re ]]; then
#	echo "Invalid output path. Usage: $0 <file_path> <out_path> <cores>"
#	exit 1
#fi

# Check if the number of cores argument is provided
#echo -n "$3" | od -c
#if [ -z "$2" ] || ! [[ $3 =~ $re ]]; then
#	echo "Invalid number of cores. Usage: $0 <file_path> <out_path> <cores>"
#	exit 1
#fi

# Print arguments
#file=$1
#outfp=$2
#cores=$3
#echo "Input: $file"
#echo "Output: $outfp"
#echo "Cores: $cores"

# Read environment variables (same as above, but additionally TRANSPOSE_FLAG)
MINKOWSKI_P=${MINKOWSKI_P:-1.5}
DISTPOW=${DISTPOW:-1.5}
LINKAGE=${LINKAGE:-wardu}
LOG2_FLAG=""
if [ "${LOG2:-1}" = "1" ]; then LOG2_FLAG="--log2"; fi
MEAN_CENTER_FLAG=""
if [ "${MEAN_CENTER:-1}" = "1" ]; then MEAN_CENTER_FLAG="--mean-center"; fi
UNIT_VARIANCE_FLAG=""
if [ "${UNIT_VARIANCE:-1}" = "1" ]; then UNIT_VARIANCE_FLAG="--unit-variance"; fi
TRANSPOSE_FLAG=""
if [ "${TRANSPOSE:-1}" = "1" ]; then TRANSPOSE_FLAG="--transpose-last"; fi

# Check arguments
if [ -z "$1" ] || [[ $1 =~ ^[0-9]+$ ]]; then
    echo "Invalid input path. Usage: $0 <file_path> <out_path> <cores>"
    exit 1
fi
if [ -z "$2" ] || [[ $2 =~ ^[0-9]+$ ]]; then
    echo "Invalid output path. Usage: $0 <file_path> <out_path> <cores>"
    exit 1
fi
if [ -z "$3" ] || ! [[ $3 =~ ^[0-9]+$ ]]; then
    echo "Invalid number of cores. Usage: $0 <file_path> <out_path> <cores>"
    exit 1
fi

file=$1
outfp=$2
cores=$3
echo "Input: $file"
echo "Output: $outfp"
echo "Cores: $cores"

/usr/bin/time \
    ./maldi-clustering-and-imaging-main/hcdist \
    ${LOG2_FLAG} ${MEAN_CENTER_FLAG} ${UNIT_VARIANCE_FLAG} ${TRANSPOSE_FLAG} \
    --tree-flip-edge --cluster-no-merge \
    --minkowski=${MINKOWSKI_P} --distpow=${DISTPOW} --${LINKAGE} --depower --threads=$cores \
    "${file}" | tee "${outfp}"
