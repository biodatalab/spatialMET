#!/usr/bin/perl -w

use File::Basename;

$imzml_filename  = shift;
$output_basename = shift;

# print usage statement
if (!defined($imzml_filename) || $imzml_filename =~ /^-/)
{
    $program_name = basename($0);

    print STDERR "Usage: $program_name image.imzML [output_filenames_prefix]\n";
    
    exit(1);
}


# sanity check filename, replace .ibd with .imzML
$filename_orig = $imzml_filename;
if ($imzml_filename =~ s/(\.(ibd|tsv|txt|zst))+\..*$/\.imzML/i)
{
    printf STDERR "# using %s as input instead of %s\n",
        $imzml_filename, $filename_orig;
}


if (!defined($output_basename))
{
    $output_basename =  $imzml_filename;
    $output_basename =~ s/\.[^.]+$//;
}

# make sure a working version of R is installed
#$module_r_str = "module add R/4.4.0";


# Rscript doesn't load methods, but R does automatically [grumble...]

# load libraries
$r_str  = "";
$r_str .= "library(methods);";
$r_str .= "library('Cardinal');";
$r_str .= "library('tidyverse');";
$r_str .= "library('data.table');";

# configure parallel processing
$r_str .= "setCardinalBPPARAM(BPPARAM=MulticoreParam(workers=32));";

# load the imzML file
$r_str .= sprintf "imzml_data=readImzML(memory=FALSE,as='MSImagingArrays',file='%s');",
                  $imzml_filename;

# pick peaks, gap-fill picked peaks with 2nd call using peaks from 1st
$r_str .= "peaks=peakProcess(imzml_data, nchunks=32);";
$r_str .= "peaks=peakProcess(imzml_data, ref=peaks, nchunks=32);";

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
#
$r_str .= "spectra=spectra(peaks);";
$r_str .= "coords=coord(peaks);";
$r_str .= q{mz_values=mz(peaks@featureData);};
$r_str .= "mat2d=as.matrix(spectra);";


# append mz_ to m/z values to prevent Excel from potentially corrupting them
# the instrument we're using only barely has 4 decimals of precision, if that
$r_str .= q{rownames(mat2d)=sprintf(\"mz_%.4f\",round(mz_values, 4));};

# format pixel col names with their x/y coordinates
$r_str .= q{colnames(mat2d)=paste0('y',coords[['y']],'x',coords[['x']]);};

# label m/z column header
$r_str .= q{mat2d=mat2d %>% as.data.frame() %>% rownames_to_column('mz_val');};

# fwrite() is orders of magnitude faster than write.csv
$r_str .= sprintf "fwrite(mat2d,nThread=32,sep='\\t',quote=F,row.names=F,'%s.tsv')",
              $output_basename;

# output the Rscript command
$rscript_str = sprintf "Rscript -e \"%s\"", $r_str;

print "$rscript_str\n";


# compress the output file
printf "zstd -f -19 -T32 \"%s.tsv\"\n", $output_basename;
printf "rm -f \"%s.tsv\"\n", $output_basename;


# print blank line(s) to better separate concatenated scripts
printf "\n";
