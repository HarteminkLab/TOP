#! /usr/bin/env Rscript

# R script to add ChIP peak labels to combined data for a TF x cell type combo
library(optparse)
library(TOP)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

# Process the command-line arguments.
option_list <- list(
  make_option('--tf_name', action='store', default=NULL, type='character',
              help='TF name.'),
  make_option('--cell_type', action='store', default=NULL, type='character',
              help='Cell type.'),
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--combined_data_dir', action='store', default=NULL, type='character',
              help='Directory for combined data for each TF x cell type combo.'),
  make_option('--chip_peak_dir', action='store', default=NULL, type='character',
              help='Directory for ChIP-seq peak bed narrowPeak files directory.'),
  make_option('--thresh_pValue', action='store', default='1e-5', type='character',
              help='FIMO p-value threshold [default: %default]'),
  make_option('--bin', action='store', default='M5', type='character',
              help='MILLIPEDE binning scheme.')
)

opt <- parse_args(OptionParser(option_list=option_list))
tf_name                <- opt$tf_name
cell_type              <- opt$cell_type
metadata_file          <- opt$metadata
combined_data_dir      <- opt$combined_data_dir
chip_peak_dir          <- opt$chip_peak_dir
thresh_pValue          <- opt$thresh_pValue
bin_method             <- opt$bin

metadata <- read.table(metadata_file, header = T, sep = '\t', stringsAsFactors = FALSE)
metadata <- metadata[metadata$tf_name == tf_name & metadata$cell_type == cell_type, ]
metadata$pwm_name <- paste(metadata$tf_name, metadata$pwm_id, thresh_pValue, sep = '_')
combined_data_file <- file.path(combined_data_dir, paste0(metadata$pwm_name, '.', metadata$cell_type, '.', bin_method, '_bins.combined_data.rds'))
chippeak_sample <- metadata$chip_peak_file

if( !file.exists(combined_data_file) || file.size(combined_data_file) == 0 ){
  cat(paste(combined_data_file, 'file does not exist or is empty. Probably no candidate sites were found. Skipped. \n'))
}else if(is.na(chippeak_sample) || chippeak_sample == ''){
  cat('No ChIP peak sample available. Skipped.')
}else{

  cat('Load combined data ... \n')
  data <- as.data.frame(readRDS(combined_data_file))
  data.gr <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)

  cat('Load ChIP peak data ... \n')
  chippeak_file <- paste0(chip_peak_dir, '/', chippeak_sample,'.bed.gz')
  if(!file.exists(chippeak_file)){
    cat('Download ChIP peak file ... \n')
    chippeak_url <- sprintf('https://www.encodeproject.org/files/%s/@@download/%s.bed.gz', chippeak_sample)
    system(paste('wget', chippeak_url, '-P', chip_peak_dir))
  }

  chip_peaks <- as.data.frame(fread(chippeak_file))
  colnames(chip_peaks) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak')
  chip_peaks.gr <- makeGRangesFromDataFrame(chip_peaks, keep.extra.columns = TRUE)

  cat('Add ChIP peak information ... \n')
  overlaps <- as.data.frame(findOverlaps(data.gr, chip_peaks.gr))
  colnames(overlaps) <- c('sites', 'chip_peaks')

  data$chip_label <- 0
  data$chip_label[overlaps$sites] <- 1

  data$chip_signalValue <- NA
  data$chip_signalValue[overlaps$sites] <- chip_peaks.gr$signalValue[overlaps$chip_peaks]

  combined_data_file <- file.path(combined_data_dir,
                                  paste0(metadata$pwm_name, '.', metadata$cell_type, '.', bin_method, '_bins.chip_peaks.combined_data.rds'))

  saveRDS(data, combined_data_file)

}


