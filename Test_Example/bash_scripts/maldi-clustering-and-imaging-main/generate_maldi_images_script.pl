#!/usr/bin/perl -w

use File::Basename;

$use_time_flag   = 1;    # run pixel tree with /usr/bin/time

$usr_bin_time_str = '';
if ($use_time_flag)
{
    $usr_bin_time_str = "/usr/bin/time"
}


$tsv_filename    = shift;
$output_basename = shift;


# print usage statement
if (!defined($tsv_filename) || $tsv_filename =~ /^-/)
{
    $program_name = basename($0);

    print STDERR "Usage: $program_name cardinal_output.tsv[.zst] [output_filenames_prefix]\n";
    
    exit(1);
}


# sanity check tsv filename, replace .imzML with .tsv
$filename_orig = $tsv_filename;
if ($tsv_filename =~ s/(\.(imzML|ibd))+/\.tsv\.zst/i)
{
    printf STDERR "# using %s as input instead of %s\n",
        $tsv_filename, $filename_orig;
}

if (!defined($output_basename))
{
    $output_basename =  $tsv_filename;
    $output_basename =~ s/\.zst+$//i;
    $output_basename =~ s/\.[^.]+$//;
    
    # strip path to original file, output everything to current path
    $output_basename = basename($output_basename);
}


$output_filtered_filename   = $output_basename . '_norm_filtered.txt';
$output_pixel_tree_filename = $output_basename . '_pixels.tree';
$output_mz_tree_filename    = $output_basename . '_mz.tree';
$output_pixel_weights_filename = $output_basename . '_pixels_weights.txt';
$output_mz_weights_filename    = $output_basename . '_mz_weights.txt';


# compressed with zstd
if ($tsv_filename =~ /\.zst$/i)
{
    $filter_str  = "# filter the input data\n";
    $filter_str .= "zstd -d -c --no-progress \\\n";
    $filter_str .= "  \"$tsv_filename\" | \\\n";
    $filter_str .= "hcdist \\\n";
    $filter_str .= "  --filter-present=0.01 --filter-unlog-sd=250 \\\n";
    $filter_str .= "  --floor-lod-ub --norm-median --output-data - | \\\n";
    $filter_str .= "hcdist \\\n";
    $filter_str .= "  --filter-log2-sd=0.5 --output-data - \\\n";
    $filter_str .= "  > \"$output_filtered_filename\"\n";
}
# regular tsv file, not compressed with zstd
else
{
    $filter_str  = "# filter the input data\n";
    $filter_str .= "hcdist \\\n";
    $filter_str .= "  --filter-present=0.01,1500 --filter-unlog-sd=250 \\\n";
    $filter_str .= "  --floor-lod-ub --norm-median --output-data \\\n";
    $filter_str .= "    \"$tsv_filename\" \\\n";
    $filter_str .= "  > \"$output_filtered_filename\"\n";
}


$mz_w_str     = "# generate m/z weights\n";
$mz_w_str    .= "hcdist \\\n";
$mz_w_str    .= "  --log2 --mean-center --unit-variance \\\n";
$mz_w_str    .= "  --tree-flip-edge --nclusters=23 --cluster-no-merge \\\n";
$mz_w_str    .= "  --minkowski=1.5 --distpow=1.5 --wardu --depower --threads=32 \\\n";
$mz_w_str    .= "  --leaf-custom \\\n";
$mz_w_str    .= "    \"$output_filtered_filename\" \\\n";
$mz_w_str    .= "  > \"$output_mz_weights_filename\"\n";


$pixel_str = "# generate the pixel tree\n";
if ($usr_bin_time_str ne '')
{
    $pixel_str .= "$usr_bin_time_str \\\n";
}
$pixel_str   .= "hcdist \\\n";
$pixel_str   .= "  --log2 --mean-center --unit-variance --transpose-last \\\n";
$pixel_str   .= "  --tree-flip-edge --nclusters=22 --cluster-no-merge \\\n";
$pixel_str   .= "  --minkowski=1.5 --distpow=1.5 --wardu --depower --threads=32 \\\n";
$pixel_str   .= "  --weight-file=\"$output_mz_weights_filename\" \\\n";
$pixel_str   .= "    \"$output_filtered_filename\" \\\n";
$pixel_str   .= "  > \"$output_pixel_tree_filename\"\n";

$pixel_w_str  = "# generate pixel weights from tree\n";
$pixel_w_str .= "hcdist --read-tree --leaf-custom \\\n";
$pixel_w_str .= "    \"$output_pixel_tree_filename\" \\\n";
$pixel_w_str .= "  > \"$output_pixel_weights_filename\"\n";

$mz_str       = "# generate the m/z tree\n";
$mz_str      .= "hcdist \\\n";
$mz_str      .= "  --log2 --mean-center --unit-variance \\\n";
$mz_str      .= "  --tree-flip-edge --nclusters=23 --cluster-no-merge \\\n";
$mz_str      .= "  --minkowski=1.5 --distpow=1.5 --wardu --depower --threads=32 \\\n";
$mz_str      .= "  --weight-file=\"$output_pixel_weights_filename\" \\\n";
$mz_str      .= "    \"$output_filtered_filename\" \\\n";
$mz_str      .= "  > \"$output_mz_tree_filename\"\n";


$maldi_str    = "# generate the MALDI images\n";
$maldi_str   .= "maldi_image_from_clusters \\\n";
$maldi_str   .= "  --nclusters-pixel=22 --nclusters-mz=23 \\\n";
$maldi_str   .= "  \"$output_filtered_filename\" \\\n";
$maldi_str   .= "  \"$output_mz_tree_filename\" \\\n";
$maldi_str   .= "  \"$output_pixel_tree_filename\" \\\n";
$maldi_str   .= "  \"$output_basename\"\n";


print $filter_str;
print "\n";
print $mz_w_str;
print "\n";
print $pixel_str;
print "\n";
print $pixel_w_str;
print "\n";
print $mz_str;
print "\n";
print $maldi_str;


# print blank lines to better separate concatenated scripts
printf "\n\n\n";
