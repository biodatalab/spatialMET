# Description

Tools for efficient tree building (hcdist) and MALDI image visualization
(maldi_image_from_clusters).  Specify hcdist --help (or any invalid option)
for a brief usage statement.  Perl scripts print usage statements when too
few arguments are specified.

<BR>



# Compiling the software

Download and install libFreeImage (https://freeimage.sourceforge.io/).
Update your LD_LIBRARY_PATH, etc. as appropriate for your compiler and
operating system, so that the C compiler can locate the .h/.a/.dll files.
hcdist does not require libFreeImage, as it handles purely text files.
<BR>

sh cchcdist.sh; sh ccmaldi.sh
<BR>

This will compile the two programs with gcc.
Edit as appropriate for other compilers.
<BR>

If you have difficulty compiling the software, stable statically-linked LINUX
binaries are located in /share/data2/welshea/bin/ on the Moffitt HPC cluster.

<BR>


# Generate pipeline scripts
```
find imzML/ -iname "*.imzML" -print0 | \
  sort -z | xargs -0 -L1 generate_cardinal_script.pl \
  > go_cardinal_all.sh

find imzML/ -iname "*.imzML" -print0 | \
  sort -z | xargs -0 -L1 generate_maldi_images_script.pl \
  > go_hcdist_all.sh
```

<BR>


# Example pipeline for MALDI visualization

## R script for Cardinal processing
```
library(methods)
library('Cardinal')
library('tidyverse')
library('data.table')

# configure parallel processing
setCardinalBPPARAM(BPPARAM=MulticoreParam(workers=32))

# load the imzML file
imzml_data=readImzML(memory=FALSE,as='MSImagingArrays',file='imzML/kras_luad_NoNormalization.imzML')

# pick peaks, gap-fill picked peaks with 2nd call using peaks from 1st
peaks=peakProcess(imzml_data, nchunks=32)
peaks=peakProcess(imzml_data, ref=peaks, nchunks=32)

# convert the S4 data structures into data structures Rscript won't puke on
#
# Rscript -e and R -e are "buggy", throw errors when R -f and R < do not
# For some reason, S4 class coersion doesn't work properly with -e ??
#  example: "no method for coercing this S4 class to a vector"
#
# We'll need to use the Cardinal mz() function to de-S4 mz first,
# otherwise, round() will error with non-numeric, and paste0() won't coerce
# into a vector.  Everything works fine in interactive R or R < script.R,
# it is only R -e and Rscript -e that break [grumble...]
spectra=spectra(peaks)
coords=coord(peaks)
mz_values=mz(peaks@featureData)
mat2d=as.matrix(spectra)

# append mz_ to m/z values to prevent Excel from potentially corrupting them
# the instrument we're using only barely has 4 decimals of precision, if that
rownames(mat2d)=sprintf(\"mz_%.4f\",round(mz_values, 4))

# format pixel col names with their x/y coordinates
colnames(mat2d)=paste0('y',coords[['y']],'x',coords[['x']])
mat2d=mat2d %>% as.data.frame() %>% rownames_to_column('mz_val')

# fwrite() is orders of magnitude faster than write.csv
fwrite(mat2d,nThread=32,sep='\t',quote=F,row.names=F,'imzML/kras_luad_NoNormalization.tsv')"
```

## compress unfiltered .tsv file, since it is generally large
```
zstd -19 -T32 \"imzML/kras_luad_NoNormalization.tsv\"
rm \"imzML/kras_luad_NoNormalization.tsv\"
```

## pre-process data file (impute, filter, normalize)
```
# filter the input data
zstd -d -c --no-progress \
  "imzML/20240130_flores_k-lung_raw.tsv.zst" | \
hcdist \
  --filter-present=0.01 --filter-unlog-sd=250 \
  --floor-lod-ub --norm-median --output-data - | \
hcdist \
  --filter-log2-sd=0.5 --output-data - \
  > "20240130_flores_k-lung_raw_norm_filtered.txt"

# generate m/z weights
hcdist \
  --log2 --mean-center --unit-variance \
  --tree-flip-edge --nclusters=23 --cluster-no-merge \
  --minkowski=1.5 --distpow=1.5 --wardu --depower --threads=48 \
  --leaf-custom \
    "20240130_flores_k-lung_raw_norm_filtered.txt" \
  > "20240130_flores_k-lung_raw_mz_weights.txt"

# generate the pixel tree
/usr/bin/time \
hcdist \
  --log2 --mean-center --unit-variance --transpose-last \
  --tree-flip-edge --nclusters=22 --cluster-no-merge \
  --minkowski=1.5 --distpow=1.5 --wardu --depower --threads=48 \
  --weight-file="20240130_flores_k-lung_raw_mz_weights.txt" \
    "20240130_flores_k-lung_raw_norm_filtered.txt" \
  > "20240130_flores_k-lung_raw_pixels.tree"

# generate pixel weights from tree
hcdist --read-tree --leaf-custom \
    "20240130_flores_k-lung_raw_pixels.tree" \
  > "20240130_flores_k-lung_raw_pixels_weights.txt"

# generate the m/z tree
hcdist \
  --log2 --mean-center --unit-variance \
  --tree-flip-edge --nclusters=23 --cluster-no-merge \
  --minkowski=1.5 --distpow=1.5 --wardu --depower --threads=48 \
  --weight-file="20240130_flores_k-lung_raw_pixels_weights.txt" \
    "20240130_flores_k-lung_raw_norm_filtered.txt" \
  > "20240130_flores_k-lung_raw_mz.tree"
```

## generate the MALDI images, heatmap, output clusters, etc..
```
maldi_image_from_clusters \
  --nclusters-pixel=22 --nclusters-mz=23 \
  "20240130_flores_k-lung_raw_norm_filtered.txt" \
  "20240130_flores_k-lung_raw_mz.tree" \
  "20240130_flores_k-lung_raw_pixels.tree" \
  "20240130_flores_k-lung_raw"
```


# m/z Peak Annotation

## Prepare HMDB/LipidMaps/DarkChem annotation resources
Install CPAN Text::Unidecode dependency: https://metacpan.org/dist/Text-Unidecode
<BR>
Download HMDB All Metabolites XML file: https://www.hmdb.ca/downloads
<BR>
Download LipidMaps SDF file: https://www.lipidmaps.org/databases/lmsd/download
<BR>
Download and install DarkChem: https://github.com/pnnl/darkchem
<BR>
<BR>
Run the following scripts to generate the input annotation tables:

```
# parse HMDB annotation, then trim it down to fewer columns with short headers
# will require ~8 gigabytes of free RAM for the two scripts piped together
parse_hmdb_xml.pl hmdb_metabolites.xml | parsed_hmdb_to_mapping_table.pl - > hmdb_mapping_table.txt

# parse LipidMaps annotation
lipidmaps_sdf_to_table.pl structures.sdf > lipidmaps_sdf_parsed.txt

# add DarkChem predictions to annotation tables
smiles_or_inchi_to_ccs.pl hmdb_mapping_table.txt       > hmdb_ccs.txt
smiles_or_inchi_to_ccs.pl lipidmaps_sdf_parsed.txt.txt > lmaps_ccs.txt

```

## Annotate m/z peaks with HMDB / LipidMaps matches

```
# convert semicolon delimited SCiLS export to tab-delimited
# SCiLS fails to escape embedded ; so we need to fix that with another script
csv2tab_not_excel.pl --semicolon 20240130_flores_k-lung_raw_mz_SCiLS.tsv | \
fix_corrupt_scils_tsv.pl - > 20240130_flores_k-lung_raw_mz_SCiLS.txt

# match HMDB and LipidMaps peaks to SCiLS and hcdist pipeline peaks
# output annotated SCiLS peaks if only SCiLS file is specified
# output annotated hcdist pipeline peaks if both are specified
#   (annotated SCiLS table will not be output in this case)
annotate_maldi.pl hmdb_ccs.txt lmaps_ccs.txt 20240130_flores_k-lung_raw_mz_SCiLS.txt 20240130_flores_k-lung_raw_mz_clusters.txt > 20240130_flores_k-lung_raw_mz_annotation.txt
```
