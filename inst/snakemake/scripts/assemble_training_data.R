#! /usr/bin/env Rscript

# Assemble training data for selected TF x cell type combos.
library(TOP)

library(optparse)

# Process the command-line arguments.
option_list <- list(
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--combined_data_dir', action='store', default=NULL, type='character',
              help='Directory of combined data for each TF x cell type combo.'),
  make_option("--thresh_pValue", action="store", default='1e-5', type='character',
              help="FIMO p-value threshold [default: %default]"),
  make_option("--bin", action="store", default='M5', type='character',
              help="MILLIPEDE binning scheme."),
  make_option("--transform", action="store", default='asinh', type='character',
              help="Transformation of ChIP counts. (default: asinh)"),
  make_option("--max_sites", action="store", default=50000, type='integer',
              help="Max number of candidate sites"),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.')
)

opt <- parse_args(OptionParser(option_list=option_list))
metadata_file          <- opt$metadata
combined_data_dir      <- opt$combined_data_dir
thresh_pValue          <- opt$thresh_pValue
bin_method             <- opt$bin
transform_chip         <- opt$transform
max_sites              <- opt$max_sites
outdir                 <- opt$outdir
outname                <- opt$outname


## Example
metadata_file          <- "/datacommons/harteminklab/kl124/TOP/data/ENCODE/metadata/processed_JASPARver1/hg38/ATAC_ChIP_JASPARver1_training_data_table.tsv"
combined_data_dir      <- "/hpc/home/kl124/work/TOP/processed_data/combined_atac_chip_data/hg38/JASPARver1"
thresh_pValue          <- "1e-5"
bin_method             <- "M5"
outdir                 <- "/hpc/home/kl124/work/TOP/processed_data/assembled_training_data/test"
outname                <- "ATAC_hg38_training_data"

##### Begins here #####

metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = FALSE)
tf_list <- sort(unique(metadata$tf_name))
celltype_list <- sort(unique(metadata$cell_type))

###################################
# Test a few TF x cell type combos
tf_list <- tf_list[1:10]
celltype_list <- celltype_list[1:3]
###################################

metadata <- metadata[which( (metadata$tf_name %in% tf_list) & (metadata$cell_type %in% celltype_list) ), ]
metadata$pwm_name <- paste(metadata$tf_name, metadata$pwm_id, thresh_pValue, sep = '_')

metadata$data_file <- file.path(combined_data_dir,
                                paste0(metadata$pwm_name, '.', metadata$cell_type, '.', bin_method, '_bins.combined_data.rds'))

tf_cell_table <- metadata[, c('tf_name', 'cell_type', 'data_file')]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.table(tf_cell_table, file.path(outdir, 'tf_cell_table.txt'), sep = '\t',
            col.names = TRUE, row.names = FALSE, quote = FALSE)

## Load data, split into partitions, and combine all tf x cell type combos
assemble_TOP_training_data(tf_cell_table_file = file.path(outdir, 'tf_cell_table.txt'),
                           training_data_dir = outdir,
                           training_data_name = 'ATAC_hg38_training_data',
                           chip_colname = 'chip',
                           training_chrs = paste0('chr', seq(1,21,2)),
                           n.partitions = 10)

cat("Assemble training data saved at:", outdir, "\n")

