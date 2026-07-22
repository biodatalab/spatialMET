##
# spatialMET (Stage 1): Processing imzML and ibd files (via Cardinal)
#

# Process user arguments
cmdargs = commandArgs(trailingOnly=TRUE)

if(any(!grepl("^\\-", cmdargs))){
  stop("Unrecognized arguments passed.")
}

cores = grep('\\-c\\=[0-9]{0,3}', cmdargs, value=T)
cores = as.numeric(gsub('\\-c\\=', '', cores))

fp = grep('\\-i\\=.+', cmdargs, value=T)
fp = gsub('\\.imzML\\s*$|\\.ibd\\s*$', '', fp)
fp = as.character(gsub('\\-i\\=', '', fp))

cat('User arguments:\n')
cat(paste0("-c: ", cores, '\n'))
cat(paste0("-i: ", fp, '\n'))

# Load libraries
library('Cardinal')
library('data.table')
library('tibble')
library('magrittr')

setCardinalBPPARAM(BPPARAM=MulticoreParam(workers=cores))

# Read imzML file 
cat('Reading imzML/ibd...\n')
start_t = Sys.time()
imzml_data=readImzML(memory=FALSE, as='MSImagingArrays', file=fp)
end_t = difftime(Sys.time(), start_t, units='min')
cat(paste0('Read-in completed in ', round(end_t, 2), ' min.\n'))

# Select peaks and align across spectra
cat('Selecting and aligning peaks...\n')
start_t = Sys.time()
peaks=peakProcess(imzml_data, nchunks=cores)
end_t = difftime(Sys.time(), start_t, units='min')
cat(paste0('Peak selection completed in ', round(end_t, 2), ' min.\n'))
start_t = Sys.time()
peaks=peakProcess(imzml_data, ref=peaks, nchunks=cores)
end_t = difftime(Sys.time(), start_t, units='min')
cat(paste0('Peak alignment completed in ', round(end_t, 2), ' min.\n'))

# Extract peak m/z values
mz_values=mz(peaks@featureData)

# Extract peak intensities
cat('Extracting intensities and coordinates...\n')
start_t = Sys.time()
spectra=spectra(peaks)
# Ensure intensities are in matrix format
mat2d=as.matrix(spectra)
# Provide row names as peak m/z values rounded to 4 digits
rownames(mat2d) = sprintf("mz_%.4f", round(mz_values, 4))

# Extract pixel coordinates
coords=coord(peaks)
coords = as.data.frame(coords)

# Add column names to intensity matrix based on coordinates
colnames(mat2d) = paste0('y',coords[['y']],'x',coords[['x']])

# Add names to pixels matching column names in intensities
coords[['pixel_name']] = paste0('y',coords[['y']],'x',coords[['x']])
coords = coords[, c('pixel_name', 'x', 'y')]

# Make rownames as column in intensity matrix for saving
mat2d = as.data.frame(mat2d) %>% 
  tibble::rownames_to_column('mz_val')

end_t = difftime(Sys.time(), start_t, units='min')
cat(paste0('Data extraction completed in ', round(end_t, 2), ' min.\n'))

# Write data to files
cat('Writing data to files...\n')
start_t = Sys.time()
fwrite(mat2d, nThread=cores, sep='\t', quote=F, row.names=F, file=paste0(fp, '_intx.tsv'))
fwrite(coords, nThread=cores, sep='\t', quote=F, row.names=F, file=paste0(fp, '_xy.tsv'))
end_t = difftime(Sys.time(), start_t, units='min')
cat(paste0('Write-in data to file completed in ', round(end_t, 2), ' min.\n'))

