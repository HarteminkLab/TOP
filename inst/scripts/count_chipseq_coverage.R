#! /usr/bin/env Rscript

# R script to count ChIP-seq read coverage for candidate sites
library(optparse)
library(TOP)

# Process the command-line arguments.
option_list <- list(
  make_option('--pwm_id', action='store', default=NULL, type='character',
              help='PWM ID.'),
  make_option('--cell_type', action='store', default=NULL, type='character',
              help='Cell type.'),
  make_option('--sites', action='store', default=NULL, type='character',
              help='Filename for candidate sites (BED format).'),
  make_option('--chrom_size', action='store', default=NULL, type='character',
              help='Filename listing the size for each chromosome.'),
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--chip_dir', action='store', default=NULL, type='character',
              help='Directory for ChIP-seq bam files'),
  make_option('--ref_size', action='store', default=2e7, type='integer',
              help='Reference library size for ChIP-seq.'),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.'),
  make_option('--path_bedtools', action='store', default='bedtools', type='character',
              help='Path to bedtools executable.')
)

opt <- parse_args(OptionParser(option_list=option_list))
pwm_id                 <- opt$pwm_id
cell_type              <- opt$cell_type
sites_file             <- opt$sites
chrom_size_file        <- opt$chrom_size
metadata_file          <- opt$metadata
chip_dir               <- opt$chip_dir
ref_size               <- opt$ref_size
outdir                 <- opt$outdir
outname                <- opt$outname
path_bedtools          <- opt$path_bedtools


# Counting ChIP-seq read coverage for candidate sites
chip_counts_file <- file.path(outdir, paste0(outname, '.normalized_chipcounts.rds'))

if( !file.exists(sites_file) ){
  cat(paste(sites_file, 'file does not exist. Probably no candidate sites were found.\n'))
  # res <- file.create(chip_counts_file) # create an empty file for snakemake downstream pipeline.
}else{
  training_metadata <- read.table(metadata_file, header = TRUE,
                                  sep = "\t", stringsAsFactors = FALSE)
  chip_counts <- count_normalize_chip(pwm_id, cell_type, training_metadata,
                                      sites_file, chrom_size_file,
                                      chip_dir, outdir,
                                      outname, path_bedtools)
  saveRDS(chip_counts, chip_counts_file)

}
