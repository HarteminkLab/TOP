#! /usr/bin/env Rscript

# Assemble training data for selected TF x cell type combos.

library(optparse)

# Process the command-line arguments.
option_list <- list(
  make_option('--tf_list', action='store', default=NULL, type='character',
              help='List of training TFs.'),
  make_option('--celltype_list', action='store', default=NULL, type='character',
              help='List of training cell types.'),
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
tf_list                <- opt$tf_list
celltype_list          <- opt$celltype_list
metadata_file          <- opt$metadata
combined_data_dir      <- opt$combined_data_dir
thresh_pValue          <- opt$thresh_pValue
bin_method             <- opt$bin
transform_chip         <- opt$transform
max_sites              <- opt$max_sites
outdir                 <- opt$outdir
outname                <- opt$outname


## Example
tf_list                <- "MITF,ZBTB7B,ZNF740,SP1,GATA4,MAFF,FOXA1,ETV1,HSF1,MXI1,ZKSCAN1,FOXJ2,GATA1"
celltype_list          <- "HepG2,MCF-7,IMR-90,A549,K562"
metadata_file          <- "/datacommons/harteminklab/kl124/TOP/data/ENCODE/metadata/processed_JASPARver1/hg38/ATAC_ChIP_JASPARver1_training_data_table.tsv"
combined_data_dir      <- "/hpc/home/kl124/work/TOP/processed_data/combined_atac_chip_data/hg38/JASPARver1"
thresh_pValue          <- "1e-5"
bin_method             <- "M5"
outdir                 <- "/hpc/home/kl124/work/TOP/processed_data/assembled_training_data/test"
outname                <- "ATAC_hg38_training_data"

##### Begins here #####

metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = FALSE)

if( !is.null(tf_list) ) {
  tf_list <- unlist(strsplit(tf_list, split = ','))
}else{
  tf_list <- sort(unique(metadata$tf_name))
}


if( !is.null(celltype_list) ) {
  celltype_list <- unlist(strsplit(celltype_list, split = ','))
}else{
  celltype_list <- sort(unique(metadata$cell_type))
}

###################################
# Test a few TF x cell type combos
tf_list <- tf_list[1:10]
celltype_list <- celltype_list[1:3]
###################################

metadata <- metadata[which( (metadata$tf_name %in% tf_list) & (metadata$cell_type %in% celltype_list) ), ]
metadata$pwm_name <- paste(metadata$tf_name, metadata$pwm_id, thresh_pValue, sep = '_')

metadata$data_file <- file.path(combined_data_dir,
                                paste0(metadata$pwm_name, '.', metadata$cell_type, '.', bin_method, '_bins.combined_data.rds'))

tf_cell_table <- metadata[, c('tf_name', 'pwm_name', 'cell_type', 'data_file')]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


## Remove missing data
missing_tf_cell_table <- NULL
rows_missing <- NULL
for(i in 1:nrow(tf_cell_table)) {
  tf_name <- as.character(tf_cell_table$tf_name[i])
  cell_type <- as.character(tf_cell_table$cell_type[i])
  data_file <- as.character(tf_cell_table$data_file[i])

  if(!file.exists(data_file) || file.size(data_file) == 0){
    cat('Warning:', tf_name, 'in', cell_type, 'data file does not exist or is empty!\n')
    cat('Check:', data_file, '\n')
    rows_missing <- c(rows_missing, i)
    missing_tf_cell_table <- rbind(missing_tf_cell_table, tf_cell_table[i,])
  }
}

if(length(rows_missing)>0){
  tf_cell_table <- tf_cell_table[-rows_missing,]
  write.table(missing_tf_cell_table, file.path(outdir, paste0(outname, '_missing.txt')),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}

tf_list <- sort(unique(tf_cell_table$tf_name))
celltype_list <- sort(unique(tf_cell_table$cell_type))
cat("TFs:", tf_list, "\n")
cat("Cell types:", celltype_list, "\n")

tf_cell_table$tf_name <- factor(tf_cell_table$tf_name, levels = tf_list)
tf_cell_table$cell_type <- factor(tf_cell_table$cell_type, levels = celltype_list)
tf_cell_table <- tf_cell_table[with(tf_cell_table, order(tf_name, cell_type)),]

fwrite(tf_cell_table, file.path(outdir, 'tf_cell_table.txt'), sep = '\t')

## Load data, split into partitions, and combine all tf x cell type combos
library(TOP)
all_training_data <- assemble_TOP_training_data(tf_cell_table = './tf_cell_table.txt',
                                                training_data_dir = './',
                                                training_data_name = 'TOP_training_data',
                                                chip_colname = 'chip',
                                                training_chrs = paste0('chr', seq(1,21,2)),
                                                n.partitions = 10)

cat("Assemble training data saved at:", outdir, "\n")

