#! /usr/bin/env Rscript

# Get DNase-seq (or ATAC-seq) count matrices for each motif
library(optparse)
library(TOP)

# Process the command-line arguments.
option_list <- list(
  make_option('--sites', action='store', default=NULL, type='character',
              help='Filename for candidate sites.'),
  make_option('--dnase_fwd', action='store', default=NULL, type='character',
              help='Filename for DNase counts in forward strand'),
  make_option('--dnase_rev', action='store', default='samtools', type='character',
              help='Filename for DNase counts in reverse strand.'),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.'),
  make_option('--path_bwtool', action='store', default='bwtool', type='character',
              help='Path to bwtool executable.')
)

opt <- parse_args(OptionParser(option_list=option_list))
file_sites           <- opt$sites
file_dnase_count_fwd <- opt$dnase_fwd
file_dnase_count_rev <- opt$dnase_rev
outdir               <- opt$outdir
outname              <- opt$outname
path_bwtool          <- opt$path_bwtool

dir.create(outdir, showWarnings = F, recursive = T)

file_dnase_matrix_fwd <- file.path(outdir, paste0(outname, '.fwd.countmatrix.txt.gz'))
file_dnase_matrix_rev <- file.path(outdir, paste0(outname, '.rev.countmatrix.txt.gz'))
file_dnase_matrix_combined <- file.path(outdir, paste0(outname, '.countmatrix.rds'))

if( !file.exists(file_sites) ){
  cat(paste(file_sites, 'file does not exist.
            Probably no candidate sites were found.\n'))
  # res <- file.create(file_dnase_matrix_combined) # create an empty file for snakemake downstream pipeline.

}else{

  # Get DNase count matrices for candidate sites
  dnase_counts <- get_dnase_sites_counts(file_sites,
                                         file_dnase_count_fwd, file_dnase_count_rev,
                                         file_dnase_matrix_fwd, file_dnase_matrix_rev,
                                         path_bwtool)

  saveRDS(dnase_counts, file_dnase_matrix_combined)

  unlink(file_dnase_matrix_fwd); unlink(file_dnase_matrix_rev)

}
