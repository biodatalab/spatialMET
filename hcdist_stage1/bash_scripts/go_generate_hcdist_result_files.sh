#!/bin/bash

##
# Generate tabular and graphical results from hcdist clustering
# Input is a text file containing a table of peak intensities
# Second argument is a string specifying the the metabolite tree.
# Third argument is a string specifying the the pixel tree.
# Fourth and Fifth arguments specify the numebr of cores to use and
# number of spatial domaind to define
#
# go_generate_hcdist_result_files.sh <intenaity_tsv> <feature_tree> <pixel_tree> <cores> <n_clusters>
#

# Check if the intesity matrix file path argument is provided
re='^[0-9]+$'
if [ -z "$1" ] || [[ $1 =~ $re ]]; then
	echo "Invalid intensity matrix file path. Usage: $0 <intx_path> <feat_path> <pix_path> <out_path> <mz_clusters> <pix_clusters>"
	exit 1
fi

# Check if the feature cluster file path argument is provided
if [ -z "$2" ] || [[ $2 =~ $re ]]; then
	echo "Invalid feature cluster file path. Usage: $0 <intx_path> <feat_path> <pix_path> <out_path> <mz_clusters> <pix_clusters>"
	exit 1
fi

# Check if the pixel cluster file path argument is provided
if [ -z "$3" ] || [[ $3 =~ $re ]]; then
	echo "Invalid pixel cluster file path. Usage: $0 <intx_path> <feat_path> <pix_path> <out_path> <mz_clusters> <pix_clusters>"
	exit 1
fi

# Check if token for result files is provided
if [ -z "$4" ] || [[ $4 =~ $re ]]; then
	echo "Missing token for result file names. Usage: $0 <intx_path> <feat_path> <pix_path> <out_token> <mz_clusters> <pix_clusters>"
	exit 1
fi

# Check if the number of predicted metabolite clusters argument is provided
if [ -z "$5" ] || ! [[ $5 =~ $re ]]; then
	echo "Invalid number of metabolite clusters. Usage: $0 <intx_path> <feat_path> <pix_path> <out_path> <mz_clusters> <pix_clusters>"
	exit 1
fi

# Check if the number of predicted pixel clusters argument is provided
if [ -z "$6" ] || ! [[ $6 =~ $re ]]; then
        echo "Invalid number of pixel clusters. Usage: $0 <intx_path> <feat_path> <pix_path> <out_path> <mz_clusters> <pix_clusters>"
        exit 1
fi

# Print arguments
intxfp=$1
featfp=$2
pixfp=$3
outfp=$4
mzclust=$5
pixclust=$6
echo "Intensities: $intxfp"
echo "Features: $featfp"
echo "Pixels: $pixfp"
echo "Output token: $outfp"
echo "Metabolite clusters: $mzclust"
echo "Predicted clusters: $pixclust"

# Generate hcdist result files
./maldi-clustering-and-imaging-main/maldi_image_from_clusters \
  --nclusters-pixel=$pixclust --nclusters-mz=$mzclust \
  "${intxfp}" \
  "${featfp}" \
  "${pixfp}" \
  "${outfp}"

