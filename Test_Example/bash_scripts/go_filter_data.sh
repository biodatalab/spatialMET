#!/bin/bash

##
# Decompress a zst file and perform filtering of the intensity data
# Input is a text file with a imzML file and number of cores to use during Cardinal 
# processing of each file, per line:
# 
# /home/user01/mass_spec_im_file_sample01.tsv.zst
# /home/user01/mass_spec_im_file_sample02.tsv.zst
# 

# Check if the file path argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <file_path>"
  exit 1
fi

# Decompress on the file with zstd and use hcdist to filter data
cat $1 | while read file; do
	outfp=$( echo $file | sed 's/\.zst$//' )
	zstd -d -c --no-progress "$file" | \
	./maldi-clustering-and-imaging-main/hcdist \
	--filter-present=0.01 \
	--filter-unlog-sd=250 \
	--floor-lod-ub --norm-median \
	--output-data - | ./maldi-clustering-and-imaging-main/hcdist \
	--filter-log2-sd=0.5 --output-data > "${outfp}_filtered.tsv"
done

