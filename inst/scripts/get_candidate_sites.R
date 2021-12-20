#! /usr/bin/env Rscript

# Get candidate TF binding sites from FIMO motif match result
library(optparse)
library(TOP)
options(scipen=999) # suppress scientific notations

# Process the command-line arguments.
option_list <- list(
  make_option("--fimo", action="store", default=NULL, type='character',
              help="Filename of FIMO result"),
  make_option("--out", action="store", default=NULL, type='character',
              help="Output filename"),
  make_option("--flank", action="store", default=100, type='integer',
              help="Size (bp) of flanking region on each side of motif [default: %default]"),
  make_option("--thresh_pValue", action="store", default=1e-5, type='double',
              help="FIMO p-value threshold [default: %default]"),
  make_option("--blacklist", action="store", default=NULL, type='character',
              help="Filename of the blacklist regions."),
  make_option("--mapability", action="store", default=NULL, type='character',
              help="Filename of the mapability reference file in bigWig format."),
  make_option("--bigWigAverageOverBed_path", action="store", default="bigWigAverageOverBed", type='character',
              help="Path to bigWigAverageOverBed executable. Only needed for computing mapability.")
  )

opt <- parse_args(OptionParser(option_list=option_list))
fimo_file                  <- opt$fimo
out_file                   <- opt$out
flank                      <- opt$flank
thresh_pValue              <- opt$thresh_pValue
blacklist_file             <- opt$blacklist
mapability_file            <- opt$mapability
bigWigAverageOverBed_path  <- opt$bigWigAverageOverBed_path

if( !file.exists(fimo_file) || file.size(fimo_file) == 0 ){
  cat(paste(fimo_file, 'file does not exist or is empty.
            Probably no motif matches were found.\n'))
  # res <- file.create(out_file) # create an empty file for snakemake downstream pipeline.

}else{

  # Get candidate sites using FIMO motif matches with flanking regions
  sites.df <- flank_fimo_sites(fimo_file, flank)

  # Select candidate sites with FIMO p-value < thresh_pValue
  sites.df <- sites.df[which(as.numeric(sites.df$p.value) < as.numeric(thresh_pValue)), ]
  cat(nrow(sites.df), 'sites with p-value <', thresh_pValue, '\n')

  # Filter candidate sites in ENCODE blacklist
  if(!is.null(blacklist_file)) {
    sites.df <- filter_blacklist(sites.df, blacklist_file)
  }

  # Compute mapability for candidate sites
  if(!is.null(mapability_file)) {
    sites.df$mapability <- compute_mapability(sites.df, mapability_file,
                                              bigWigAverageOverBed_path)
  }

  if(!dir.exists(dirname(out_file))){
    dir.create(dirname(out_file), showWarnings = F, recursive = T)
  }

  colnames(sites.df)[1] <- paste0('#', colnames(sites.df)[1])

  write.table(sites.df, out_file, sep = '\t',
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  cat(nrow(sites.df), 'candidate sites written in', out_file, '\n')
}

