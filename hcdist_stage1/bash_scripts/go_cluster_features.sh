#!/bin/bash

##
# Cluster features (i.e., metabolites) from intensity data using hcdist
# Input is a text file containing a table of peak intensities
# Second argument is a string specifying the name of output file containing 
# the metabolite tree.
# Third argument specifies the number of cores to use.
#
# go_cluster_features.sh /home/user01/intensity_table.tsv_filtered.tsv /home/user01/intensity_table_feature_tree.tre 32
# 

# Check if the file path argument is provided
re='^[0-9]+$'
if [ -z "$1" ] || [[ $1 =~ $re ]]; then
	echo "Invalid input path. Usage: $0 <file_path> <out_path> <cores>"
	exit 1
fi

# Check if the output file path s provided
if [ -z "$2" ] || [[ $2 =~ $re ]]; then
	echo "Invalid output path. Usage: $0 <file_path> <out_path> <cores>"
	exit 1
fi

# Check if the number of cores argument is provided
#echo -n "$3" | od -c
if [ -z "$3" ] || ! [[ $3 =~ $re ]]; then
	echo "Invalid number of cores. Usage: $0 <file_path> <out_path> <cores>"
	exit 1
fi

# Print arguments
file=$1
outfp=$2
cores=$3
echo "Input: $file"
echo "Output: $outfp"
echo "Cores: $cores"

# Cluster features using hcdist
/usr/bin/time \
	./maldi-clustering-and-imaging-main/hcdist \
	--log2 --mean-center --unit-variance \
	--tree-flip-edge --cluster-no-merge \
	--minkowski=1.5 --distpow=1.5 --wardu --depower --threads=$cores \
	"${file}" | tee "${outfp}"

