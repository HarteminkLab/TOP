#! /usr/bin/env Rscript

# R script to count DNase (or ATAC) cuts around candidate sites.
library(optparse)
library(TOP)

# Process the command-line arguments.
option_list <- list(
  make_option("--bam", action="store", default=NULL, type='character',
              help="BAM filename of DNase-seq alignments."),
  make_option("--chrom_size", action="store", default=NULL, type='character',
              help="Filename listing the size for each chromosome."),
  make_option("--out_dir", action="store", default=NULL, type='character',
              help="Output directory."),
  make_option("--out_prefix", action="store", default=NULL, type='character',
              help="Output filename prefix."),
  make_option("--bedtools_path", action="store", default='bedtools', type='character',
              help="Path to bedtools executable."),
  make_option("--bedGraphToBigWig_path", action="store", default='bedGraphToBigWig', type='character',
              help="Path to UCSC bedGraphToBigWig executable."),
  make_option("--bedSort_path", action="store", default='bedSort', type='character',
              help="Path to UCSC bedSort executable.")
)

opt <- parse_args(OptionParser(option_list=option_list))
bam_file               <- opt$bam
chrom_size_file       <-  opt$chrom_size
out_dir                <- opt$out_dir
out_prefix             <- opt$out_prefix
bedtools_path          <- opt$bedtools_path
bedGraphToBigWig_path  <- opt$bedGraphToBigWig_path
bedSort_path           <- opt$bedSort_path


# Count DNase-seq cuts (5' end) along the genome, and save in Bigwig format.
count_dnase_genome_cuts(bam_file,
                        chrom_size_file,
                        out_dir,
                        out_prefix,
                        bedtools_path,
                        bedGraphToBigWig_path)

