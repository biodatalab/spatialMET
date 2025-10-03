#!/bin/bash

##
# Serially process files (imzML+ibd) to extract data and coordinates using Cardinal
# Compresses the resulting intensity matrix (text file)
# Input is a text file with a imzML file and number of cores to use during Cardinal 
# processing of each file, per line:
# 
# /home/user01/mass_spec_im_file_sample01.imZML 32
# /home/user01/mass_spec_im_file_sample02.imZML 32
#
# HPC users: Probably need to 'module load R' before running the script
# 

# Check if the file path argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <file_path> <cores>"
  exit 1
fi

cat $1 | while read line; do
	cores=$( echo $line | grep -oE "[0-9]+$" )
	file=$( echo $line | sed -E 's/[0-9]+$//' )
	outfp=$( echo $line | sed -E 's/[0-9]+$//' | sed 's/\.imzML\s$//' )
	
	echo "Bash args:"
	echo ${cores}
	echo ${file}
	echo ${outfp}
	
	Rscript R_templates/spatialMET_stage_1_imzml_ibd_processing.R -c="${cores}" -i="${file}" -o="${outfp}"
	zstd -f -19 -T${cores} "${outfp}_intx.tsv"
	rm -f "${outfp}_intx.tsv"
done

