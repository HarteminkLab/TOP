
#' @title Sort and index the BAM file, and then retrieve and print stats in the index file.
#'
#' @param bam_file Input BAM file.
#' @param dir_output Output directory.
#' @param sort logical. If TRUE, sort (and index) the BAM file.
#' @param stats logical. If TRUE, retrieve and print stats in the index file.
#' @param path_samtools Path to samtools executable.
#'
#' @export
bam_sort_index_stats <- function(bam_file, dir_output=NA,
                                 sort=TRUE, stats=TRUE,
                                 path_samtools='samtools') {


  if ( system( paste(path_samtools,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( "ERROR: samtools could not be executed, set path_samtools\n" , sep='', file=stderr() )
    cleanup()
    q()
  }

  if( is.na(dir_output) ) {
    dir_output <- dirname(bam_file)
  } else {
    dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)
  }

  bam_prefix <- gsub('.bam$', '', basename(bam_file))
  bam_file_sorted <- paste0(dir_output, '/', bam_prefix, '.bam')

  if( sort == TRUE ) {
    # sort and index the bam file
    cat('Sort and index the bam file...\n')
    bam_file_tmp <- paste0(dir_output, '/', bam_prefix, '.tmp.bam')

    cmd <- paste(path_samtools, 'sort', bam_file, '-o', bam_file_tmp)
    system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

    file.rename(from = bam_file_tmp, to = bam_file_sorted)

    cmd <- paste(path_samtools, 'index', bam_file_sorted)
    system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
  }

  if( stats == TRUE ){
    # Get index statistics (including the number of mapped reads in the third column)
    cat('Retrieve and print stats in the index file. Result saved to: ', bam_idxstats, '\n', sep='')
    bam_idxstats <- paste0(dir_output, '/', bam_prefix, '_bam_idxstats.txt')
    cmd <- paste(path_samtools, 'idxstats', bam_file_sorted, '>', bam_idxstats)
    system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )
  }

}


#' @title Count the coverage of DNase-seq cuts for a given genome.
#'
#' @param bam_file Input BAM file.
#' @param file_out Output file name.
#' @param strand strand: '+' or '-'.
#' @param file_genome_sizes File name of genome sizes by chromosomes.
#' @param path_samtools Path to samtools executable.
#' @param path_bedtools Path to bedtools executable.
#' @param path_bedGraphToBigWig Path to bedGraphToBigWig executable.
#' @param output_format Output format: 'Bigwig' or 'bedGraph'.
#'
#' @export
count_genome_coverage <- function(bam_file,
                                  file_out=NA,
                                  strand='+',
                                  file_genome_sizes=NA,
                                  path_samtools='samtools',
                                  path_bedtools='bedtools',
                                  path_bedGraphToBigWig='bedGraphToBigWig',
                                  output_format='Bigwig') {

  if ( system( paste(path_samtools,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( "ERROR: samtools could not be executed, set path_samtools\n" , sep='', file=stderr() )
    cleanup()
    q()
  }

  if ( system( paste(path_bedtools,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( "ERROR: bedtools could not be executed, set path_bedtools\n" , sep='', file=stderr() )
    cleanup()
    q()
  }

  bam_prefix <- gsub('.bam$', '', basename(bam_file))

  if( is.na(file_out) ) {
    file_out <- paste0(dirname(bam_file), '/', bam_prefix, '.5p.tagcount')
  }

  if( sort == TRUE ) {
    bam_sort_index_stats(bam_file, dir_output=dirname(bam_file), sort=TRUE, stats=FALSE, path_samtools)
  }

  cat('Compute genome coverage for', bam_file, 'of', strand, 'strand and output in bedGraph format...\n', sep='')

  if( output_format != 'bedGraph' ) {
    file_bg <- paste0(tools::file_path_sans_ext(file_out), '.bedGraph')
  }else{
    file_bg <- file_out
  }

  cmd <- paste(path_bedtools, 'genomecov -bg -5',
               '-strand', strand, '-ibam', bam_file, '-g', file_genome_sizes,
               '>', file_bg)
  system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

  if( output_format == 'Bigwig' ){
    cat('Convert bedGraph to Bigwig format...\n', sep='')

    if ( system( paste(path_bedGraphToBigWig,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
      cat( "ERROR: bedGraphToBigWig could not be executed, set path_bedGraphToBigWig\n" , sep='', file=stderr() )
      cleanup()
      q()
    }

    cmd <- paste(path_bedGraphToBigWig, file_bg, file_genome_sizes, file_out)
    system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

    file.remove(file_bg)
  }

}

#' @title Match DNase tagcount matrices for candidate sites.
#'
#' @param file_sites filename for candidate sites
#' @param file_dnase_tagcount_fwd filename for DNase tagcount in forward strand
#' @param file_dnase_tagcount_rev filename for DNase tagcount in reverse strand
#' @param file_dnase_matrix_fwd filename for DNase tagcount matrix in forward strand
#' @param file_dnase_matrix_rev filename for DNase tagcount matrix in reverse strand
#' @param path_bwtool Path to bwtool executable.
#'
#' @export
match_tagcount_sites <- function(file_sites,
                                 file_dnase_tagcount_fwd, file_dnase_tagcount_rev,
                                 file_dnase_matrix_fwd, file_dnase_matrix_rev,
                                 path_bwtool='bwtool') {

  cat('Match matrix of', file_dnase_matrix_fwd, '... \n')
  cmd <- paste('cut -f 1-4', file_sites, '|', path_bwtool, 'extract bed stdin',
               file_dnase_tagcount_fwd, file_dnase_matrix_fwd, '-fill=0 -decimals=0 -tabs')
  system(cmd)

  cat('Match matrix of', file_dnase_matrix_rev, '... \n')
  cmd <- paste('cut -f 1-4', file_sites, '|', path_bwtool, 'extract bed stdin',
               file_dnase_tagcount_rev, file_dnase_matrix_rev, '-fill=0 -decimals=0 -tabs')
  system(cmd)

  # Flip the tagcounts generated from bwtool for motifs on the reverse strand
  rev_tagcount_bwtool(file_sites, file_dnase_matrix_fwd, file_dnase_matrix_rev)
}

#' @title Flip the tagcounts generated from bwtool for motifs on the reverse (minus) strand
#'
#' @param file_sites filename for candidate sites
#' @param file_dnase_matrix_fwd filename for DNase tagcount matrix in forward strand
#' @param file_dnase_matrix_rev filename for DNase tagcount matrix in reverse strand
rev_tagcount_bwtool <- function(file_sites, file_dnase_matrix_fwd, file_dnase_matrix_rev) {

  cat("Flipping tagcounts for motifs on the reverse strand ... \n")

  fwd_count.df <- read.table(file_dnase_matrix_fwd)
  rev_count.df <- read.table(file_dnase_matrix_rev)
  sites.df <- read.table(file_sites)

  if(sum(sites.df[,2] != fwd_count.df[,2]) != 0) {
    stop("Sites do not match!")

  } else {

    # Extract the count values
    fwd_count.m <- fwd_count.df[, 6:ncol(fwd_count.df)]
    rev_count.m <- rev_count.df[, 6:ncol(rev_count.df)]

    fwd_output <- fwd_count.m
    rev_output <- rev_count.m

    # For motifs match to the minus strand, flip the fwd and rev tagcounts, and reverse the counts
    idx_minusStrand <- which(sites.df[,6] == "-")

    fwd_output[idx_minusStrand, ] <- t(apply(rev_count.m[idx_minusStrand, ], 1, rev))
    rev_output[idx_minusStrand, ] <- t(apply(fwd_count.m[idx_minusStrand, ], 1, rev))

    # Add the sites info to the first few columns
    fwd_count.df <- cbind(sites.df[,c(1:3,6)], fwd_output)
    rev_count.df <- cbind(sites.df[,c(1:3,6)], rev_output)

    write.table(fwd_count.df, gzfile(file_dnase_matrix_fwd), sep = " ", quote = F, row.names = F, col.names = F)
    write.table(rev_count.df, gzfile(file_dnase_matrix_rev), sep = " ", quote = F, row.names = F, col.names = F)

  }
}
