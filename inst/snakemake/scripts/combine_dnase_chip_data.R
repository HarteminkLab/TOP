#! /usr/bin/env Rscript

# R script to combine DNase and ChIP data for a TF x cell type combo
library(optparse)
library(TOP)

# Process the command-line arguments.
option_list <- list(
  make_option('--tf_name', action='store', default=NULL, type='character',
              help='TF name.'),
  make_option('--cell_type', action='store', default=NULL, type='character',
              help='Cell type.'),
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option("--thresh_pValue", action="store", default='1e-5', type='character',
              help="FIMO p-value threshold [default: %default]"),
  make_option('--chip_counts', action='store', default=NULL, type='character',
              help='Filename of ChIP counts'),
  make_option('--dnase_dir', action='store', default=NULL, type='character',
              help='Directory DNase counts. '),
  make_option('--dnase_idxstats_dir', action='store', default=NULL, type='character',
              help='Directory for DNase idxstats files.'),
  make_option('--dnase_ref_size', action='store', default=1e8, type='integer',
              help='Reference library size for DNase-seq.'),
  make_option("--bin", action="store", default='M5', type='character',
              help="MILLIPEDE binning scheme."),
  make_option("--transform", action="store", default='asinh', type='character',
              help="Transformation of binned DNase counts. (asinh or log2)"),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.')
)

opt <- parse_args(OptionParser(option_list=option_list))
tf_name                <- opt$tf_name
cell_type              <- opt$cell_type
metadata_file          <- opt$metadata
thresh_pValue          <- opt$thresh_pValue
chip_counts_file       <- opt$chip_counts
dnase_countmatrix_dir  <- opt$dnase_dir
dnase_idxstats_dir     <- opt$dnase_idxstats_dir
dnase_ref_size         <- opt$dnase_ref_size
bin_method             <- opt$bin
transform              <- opt$transform
outdir                 <- opt$outdir
outname                <- opt$outname


if( !file.exists(chip_counts_file) || file.size(chip_counts_file) == 0 ){
  cat(paste(chip_counts_file, 'file does not exist or is empty. Probably no candidate sites were found.\n'))
} else {

  # Select metadata for the TF x cell type combo
  training_metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = FALSE)
  selected_metadata <- training_metadata[which(training_metadata$tf_name == tf_name & training_metadata$cell_type == cell_type), ]
  pwm_id <- selected_metadata$pwm_id
  name_acc <- grep("acc_file|dnase_file|atac_file", colnames(selected_metadata), value = T, ignore.case = T)
  if(length(name_acc) > 1){
    stop('More than 1 column of file names of DNase (or ATAC)! Please check the metadata!')
  }

  dnase_samples <- unlist(strsplit(selected_metadata[, name_acc],';'))

  dir.create(outdir, showWarnings = F, recursive = T)

  if ( length(dnase_samples) == 1 ) {

    cat('Combining DNase and ChIP data of', tf_name, 'in', cell_type, '.\n')

    dnase_counts_file <- paste0(dnase_countmatrix_dir,'/', dnase_samples, '/',
                                pwm_id, '_', thresh_pValue, '.normalized_countmatrix.rds')
    combined_data <- combine_dnase_bins_chip_data(dnase_counts_file, chip_counts_file, bin_method)

    combined_data_file <- file.path(outdir, paste0(outname, '.', bin_method, '_bins.combined_data.rds'))
    saveRDS(combined_data, combined_data_file)

  } else if ( length(dnase_samples) > 1 ) {
    cat('Combining DNase and ChIP data of', tf_name, 'in', cell_type, 'with', length(dnase_samples), 'replicates...\n')

    # Combine data for DNase replicates separately
    for (i in 1:length(dnase_samples)) {
      cat('Combining DNase and ChIP data for Rep', i, ':\t', dnase_samples[i], '...\n')
      dnase_counts_file <- paste0(dnase_countmatrix_dir,'/', dnase_samples[i],
                                  '/', pwm_id, '_', thresh_pValue, '.normalized_countmatrix.rds')
      combined_data <- combine_dnase_bins_chip_data(dnase_counts_file, chip_counts_file, bin_method)

      combined_data_file <- file.path(outdir, paste0(outname, '.', bin_method, '_bins.dnase_rep', i, '.combined_data.rds'))
      saveRDS(combined_data, combined_data_file)
    }

    # Combine data by first merging DNase replicates
    cat('Merge DNase samples:', dnase_samples, '...\n')
    dnase_counts_files <- paste0(dnase_countmatrix_dir, '/', dnase_samples, '/', pwm_id, '_', thresh_pValue, '.countmatrix.rds')
    dnase_idxstats_files <- paste0(dnase_idxstats_dir, '/', dnase_samples, '.bam.idxstats.txt')

    # merging DNase replicates
    merged_dnase_data <- merge_normalize_dnase(dnase_counts_files, dnase_idxstats_files, dnase_ref_size)
    dir.create(paste0(dnase_countmatrix_dir,'/merged_reps/'), showWarnings = FALSE, recursive = TRUE)
    merged_dnase_counts_file <- paste0(dnase_countmatrix_dir,'/merged_reps/',
                                       pwm_id, '_', thresh_pValue, '.', cell_type, '.merged_normalized_countmatrix.rds')
    saveRDS(merged_dnase_data, merged_dnase_counts_file)

    combined_data <- combine_dnase_bins_chip_data(merged_dnase_counts_file, chip_counts_file, bin_method)

    combined_data_file <- file.path(outdir, paste0(outname, '.', bin_method, '_bins.combined_data.rds'))
    saveRDS(combined_data, combined_data_file)

  }

}

