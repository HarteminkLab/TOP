
#' @title Counts DNase-seq or ATAC-seq cuts along the genome
#' @description Counts genomic cuts (5' end) from DNase-seq or
#' ATAC-seq BAM alignment files using \code{bedtools}
#' For ATAC-seq, when \code{shift_ATAC = TRUE}, shifts reads
#' so as to address offsets and align the signal across strands.
#' @param bam_file Sorted BAM file.
#' @param chrom_size_file Chromosome size file.
#' @param data_type Data type. Options: \sQuote{DNase} or \sQuote{ATAC}.
#' @param shift_ATAC Logical. When \code{shift_ATAC=TRUE} (and \code{data_type='ATAC'}),
#' shifts reads according to \code{shift_ATAC_bases}.
#' @param shift_ATAC_bases Number of bases to shift on + and - strands.
#' Default: shifts reads on + strand by 4 bp and reads on - strand by -4 bp.
#' @param outdir Output directory (default: use the directory of \code{bam_file}).
#' @param outname Output prefix (default: use the prefix of \code{bam_file}).
#' @param bedtools_path Path to \code{bedtools} executable.
#' @param bedGraphToBigWig_path Path to UCSC \code{bedGraphToBigWig} executable.
#' @export
#' @examples
#' \dontrun{
#' # ATAC-seq data
#' count_genome_cuts(bam_file='K562.ATAC.bam',
#'                   chrom_size_file='hg38.chrom.sizes',
#'                   data_type='ATAC',
#'                   shift_ATAC=TRUE,
#'                   outdir='processed_data',
#'                   outname='K562.ATAC')
#'
#' # DNase-seq data
#' count_genome_cuts(bam_file='K562.DNase.bam',
#'                   chrom_size_file='hg38.chrom.sizes',
#'                   data_type='DNase',
#'                   outdir='processed_data',
#'                   outname='K562.DNase')
#'}
count_genome_cuts <- function(bam_file,
                              chrom_size_file,
                              data_type=c('DNase', 'ATAC'),
                              shift_ATAC=TRUE,
                              shift_ATAC_bases=c(4L,-4L),
                              outdir=dirname(bam_file),
                              outname,
                              bedtools_path='bedtools',
                              bedGraphToBigWig_path='bedGraphToBigWig'){

  # Checking input arguments
  if ( Sys.which(bedtools_path) == '')
    stop( 'bedtools could not be executed. Please install bedtools and set bedtools_path.' )

  # Checking input arguments
  if ( Sys.which(bedGraphToBigWig_path) == '' )
    stop( 'bedGraphToBigWig could not be executed. Please install bedGraphToBigWig and set bedGraphToBigWig_path.' )

  data_type <- match.arg(data_type)

  if(missing(outdir))
    outdir <- dirname(bam_file)

  if(!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  if(missing(outname))
    outname <- tools::file_path_sans_ext(basename(bam_file))

  for (strand in c('+', '-')) {
    if (strand == '+') {
      bedgraph_file <- file.path(outdir, paste0(outname, '.fwd.genomecounts.bedGraph'))
    } else {
      bedgraph_file <- file.path(outdir, paste0(outname, '.rev.genomecounts.bedGraph'))
    }

    # Count 5' end cuts and save as a .bedGraph file
    cat('Counting genome cuts for', bam_file, 'on', strand, 'strand...\n')
    cmd <- paste(bedtools_path, 'genomecov -bg -5', '-strand', strand,
                 '-ibam', bam_file, '>', bedgraph_file)
    if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

    cat('Sorting counted genome cuts on', strand, 'strand...\n')
    genome_counts <- data.table::fread(bedgraph_file)
    genome_counts <- genome_counts[with(genome_counts, order(V1, V2)), ]

    # Shift ATAC-seq reads to get the centers of Tn5 binding positions
    if (data_type == 'ATAC' && shift_ATAC == TRUE) {
      if (strand == '+') {
        cat(sprintf('Shifting ATAC-seq reads by %d bp ...\n', shift_ATAC_bases[1]))
        genome_counts[,c(2:3)] <- genome_counts[,c(2:3)] + shift_ATAC_bases[1]
      } else if (strand == '-'){
        cat(sprintf('Shifting ATAC-seq reads by %d bp ...\n', shift_ATAC_bases[2]))
        genome_counts[,c(2:3)] <- genome_counts[,c(2:3)] + shift_ATAC_bases[2]
      }
    }

    data.table::fwrite(genome_counts, bedgraph_file,
                       sep='\t', col.names = FALSE, scipen = 999)

    # save to bigwig file
    bedGraphToBigWig(bedgraph_file, chrom_size_file, bedGraphToBigWig_path)
    unlink(bedgraph_file)
  }

}

#' @title Extracts count matrices around candidate binding sites
#' @description Extracts counts around candidate binding sites on both strands
#' from the genome counts data
#' (BigWig files generated using \code{count_genome_cuts()}).
#' It utilizes the \code{extract bed} function from the \code{bwtool} software
#' to extract the read counts,
#' then combines the counts into one matrix, with the first half of the columns
#' representing the read counts on the forward strand,
#' and the second half of the columns representing the read counts
#' on the reverse strand.
#' @param sites A data frame containing the candidate sites.
#' @param genomecount_dir Directory for genome counts,
#' the same as \code{outdir} in \code{count_genome_cuts()}.
#' @param genomecount_name File prefix for genome counts,
#' the same as \code{outname} in \code{count_genome_cuts()}.
#' @param tmpdir Temporary directory to save intermediate files.
#' @param bedGraphToBigWig_path Path to UCSC \code{bedGraphToBigWig} executable.
#' @param bwtool_path Path to \code{bwtool} executable.
#' @return A count matrix. The first half of the columns
#' are the read counts on the forward strand, and the second half of the
#' columns are the read counts on the reverse strand.
#' @export
#' @examples
#' \dontrun{
#' # Extracts ATAC-seq count matrices around candidate sites
#' count_matrix <- get_sites_counts(sites,
#'                                  genomecount_dir='processed_data',
#'                                  genomecount_name='K562.ATAC')
#' }
get_sites_counts <- function(sites,
                             genomecount_dir,
                             genomecount_name,
                             tmpdir=genomecount_dir,
                             bedGraphToBigWig_path='bedGraphToBigWig',
                             bwtool_path='bwtool') {

  # Checking input arguments
  if ( Sys.which(bedGraphToBigWig_path) == '' )
    stop( 'bedGraphToBigWig could not be executed. Please install bedGraphToBigWig and set bedGraphToBigWig_path.' )

  if ( Sys.which(bwtool_path) == '' ) {
    stop( 'bwtool could not be executed. Please install bwtool and set bwtool_path.' )
  }

  genome_fwd_count_file <- file.path(genomecount_dir, paste0(genomecount_name, '.fwd.genomecounts.bw'))
  genome_rev_count_file <- file.path(genomecount_dir, paste0(genomecount_name, '.rev.genomecounts.bw'))

  if( !(file.exists(genome_fwd_count_file) & file.exists(genome_rev_count_file)) ){
    stop( 'Genome count files cannot be located. Please run count_genome_cuts() first.' )
  }

  cat('Extract counts around candidate sites ... \n')
  sites_file <- tempfile(pattern = 'sites.', tmpdir = tmpdir, fileext = '.txt')
  data.table::fwrite(sites[,1:4], sites_file, sep = '\t', col.names = FALSE, scipen = 999)

  fwd_matrix_file <- tempfile(pattern = 'sitescounts.', tmpdir = tmpdir, fileext = '.fwd.matrix')
  rev_matrix_file <- tempfile(pattern = 'sitescounts.', tmpdir = tmpdir, fileext = '.rev.matrix')

  cmd <- paste(bwtool_path, 'extract bed', sites_file,
               genome_fwd_count_file, fwd_matrix_file, '-fill=0 -decimals=0 -tabs')
  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

  cmd <- paste(bwtool_path, 'extract bed', sites_file,
               genome_rev_count_file, rev_matrix_file, '-fill=0 -decimals=0 -tabs')
  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

  # Flips the counts generated from bwtool for motifs on the reverse strand and combine counts on both strands
  sites_counts.mat <- flip_neg_strand_counts(sites, fwd_matrix_file, rev_matrix_file)

  unlink(c(sites_file, fwd_matrix_file, rev_matrix_file))
  return(sites_counts.mat)
}

# Flips the counts for motif matches on the - strand
flip_neg_strand_counts <- function(sites,
                                   fwd_matrix_file,
                                   rev_matrix_file,
                                   write_updated_matrix_files = FALSE) {

  fwd_count <- as.data.frame(data.table::fread(fwd_matrix_file))
  rev_count <- as.data.frame(data.table::fread(rev_matrix_file))

  if(sum(sites[,2] != fwd_count[,2]) != 0) {
    stop('Sites do not match!')
  }

  # Extracts the count values
  fwd_count.m <- as.matrix(fwd_count[, -c(1:5)])
  colnames(fwd_count.m) <- paste0('F', c(1:ncol(fwd_count.m)))
  rev_count.m <- as.matrix(rev_count[, -c(1:5)])
  colnames(rev_count.m) <- paste0('R', c(1:ncol(rev_count.m)))

  # For motifs match to the - strand, flips the fwd and rev counts, and reverse the counts
  sites_counts.l <- list(fwd = fwd_count.m, rev = rev_count.m)
  neg_strand <- which(sites$strand == '-')
  sites_counts.l$fwd[neg_strand, ] <- t(apply(rev_count.m[neg_strand, ], 1, rev))
  sites_counts.l$rev[neg_strand, ] <- t(apply(fwd_count.m[neg_strand, ], 1, rev))

  sites_counts.mat <- as.matrix(cbind(sites_counts.l$fwd, sites_counts.l$rev))
  rownames(sites_counts.mat) <- sites$name

  if(write_updated_matrix_files){
    # writes the updated the matrix counts files
    fwd_count <- cbind(sites[,c(1:3,6)], sites_counts.l$fwd)
    rev_count <- cbind(sites[,c(1:3,6)], sites_counts.l$rev)
    data.table::fwrite(fwd_count, fwd_matrix_file, sep = ' ', scipen = 999)
    data.table::fwrite(rev_count, rev_matrix_file, sep = ' ', scipen = 999)
  }

  return(sites_counts.mat)
}

#' @title Performs \code{MILLIPEDE} binning on count matrix
#' @description Performs binning using different \code{MILLIPEDE} binning schemes
#' (M5, M24, M12, M3, M2, M1) on the input count matrix.
#' @param counts DNase-seq or ATAC-seq read counts matrix,
#' rows are candidate sites,
#' columns are DNase or ATAC counts with 100bp flanks around motifs
#' on the forward and reverse strands.
#' @param bin_method \code{MILLIPEDE} binning scheme. Options:
#' \sQuote{M5} (default), \sQuote{M24}, \sQuote{M12}, \sQuote{M3}, \sQuote{M2},
#' and \sQuote{M1}.
#' @param combine_strands Method to combine counts on both strands from M24 bins
#' to M12 bins: \sQuote{vertical} (combine counts from both strands vertically)
#' or \sQuote{motif} (combine counts from both strands with respect to
#' motif match orientation).
#' @return A list containing binning results (data frames) using different
#' binning schemes.
#' @export
#' @examples
#' \dontrun{
#' # Performs MILLIPEDE binning with different binning schemes.
#'
#' # M5 binning
#' M5_bins <- millipede_binning(counts, bin_method = 'M5')
#' }
millipede_binning <- function(counts,
                              bin_method = c('M5','M24','M12','M3','M2','M1'),
                              combine_strands=c('vertical', 'motif')) {

  bin_method <- match.arg(bin_method)
  combine_strands <- match.arg(combine_strands)

  counts <- as.matrix(counts)

  # bin locations
  flank_LF1 <- 81:100
  flank_LF2 <- flank_LF1 - 20
  flank_LF3 <- flank_LF2 - 20
  flank_LF4 <- flank_LF3 - 20
  flank_LF5 <- flank_LF4 - 20
  flank_RF1 <- (ncol(counts)/2-99): (ncol(counts)/2-80)
  flank_RF2 <- flank_RF1 + 20
  flank_RF3 <- flank_RF2 + 20
  flank_RF4 <- flank_RF3 + 20
  flank_RF5 <- flank_RF4 + 20

  flank_LR1 <- ncol(counts)/2 + flank_LF1
  flank_LR2 <- ncol(counts)/2 + flank_LF2
  flank_LR3 <- ncol(counts)/2 + flank_LF3
  flank_LR4 <- ncol(counts)/2 + flank_LF4
  flank_LR5 <- ncol(counts)/2 + flank_LF5
  flank_RR1 <- ncol(counts)/2 + flank_RF1
  flank_RR2 <- ncol(counts)/2 + flank_RF2
  flank_RR3 <- ncol(counts)/2 + flank_RF3
  flank_RR4 <- ncol(counts)/2 + flank_RF4
  flank_RR5 <- ncol(counts)/2 + flank_RF5

  motif_LF <- 101: round(ncol(counts)/4)
  motif_RF <- (round(ncol(counts)/4)+1) : (ncol(counts)/2-100)
  motif_LR <- ncol(counts)/2 + motif_LF
  motif_RR <- ncol(counts)/2 + motif_RF

  # M24 bins
  M24 <- data.frame(
    flank_LF5 <- rowSums(counts[,flank_LF5]),
    flank_LF4 <- rowSums(counts[,flank_LF4]),
    flank_LF3 <- rowSums(counts[,flank_LF3]),
    flank_LF2 <- rowSums(counts[,flank_LF2]),
    flank_LF1 <- rowSums(counts[,flank_LF1]),

    motif_LF <- rowSums(counts[,motif_LF]),
    motif_RF <- rowSums(counts[,motif_RF]),

    flank_RF1 <- rowSums(counts[,flank_RF1]),
    flank_RF2 <- rowSums(counts[,flank_RF2]),
    flank_RF3 <- rowSums(counts[,flank_RF3]),
    flank_RF4 <- rowSums(counts[,flank_RF4]),
    flank_RF5 <- rowSums(counts[,flank_RF5]),

    flank_LR5 <- rowSums(counts[,flank_LR5]),
    flank_LR4 <- rowSums(counts[,flank_LR4]),
    flank_LR3 <- rowSums(counts[,flank_LR3]),
    flank_LR2 <- rowSums(counts[,flank_LR2]),
    flank_LR1 <- rowSums(counts[,flank_LR1]),

    motif_LR <- rowSums(counts[,motif_LR]),
    motif_RR <- rowSums(counts[,motif_RR]),

    flank_RR1 <- rowSums(counts[,flank_RR1]),
    flank_RR2 <- rowSums(counts[,flank_RR2]),
    flank_RR3 <- rowSums(counts[,flank_RR3]),
    flank_RR4 <- rowSums(counts[,flank_RR4]),
    flank_RR5 <- rowSums(counts[,flank_RR5])
  )

  # M12 bins
  if(combine_strands == 'vertical'){
    M12 <- M24[,c(1:12)] + M24[,c(13:24)]
  }else if (combine_strands == 'motif') {
    M12 <- M24[,c(1:12)] + M24[,c(24:13)]
  }
  colnames(M12) <- c('L5', 'L4', 'L3', 'L2', 'L1',
                     'motif_L', 'motif_R',
                     'R1', 'R2', 'R3', 'R4', 'R5')

  # M5 bins
  left2_bin <- rowSums(M12[,c(1,2)])
  left1_bin <- rowSums(M12[,c(3,4,5)])
  motif_bin <- rowSums(M12[,c(6,7)])
  right1_bin <- rowSums(M12[,c(8,9,10)])
  right2_bin <- rowSums(M12[,c(11,12)])

  if (bin_method == 'M5') {
    bins <- data.frame(left2_bin, left1_bin, motif_bin, right1_bin, right2_bin)
  }else if (bin_method == 'M3'){
    bins <- data.frame(left1_bin, motif_bin, right1_bin)
  }else if (bin_method == 'M2'){
    bins <- data.frame(left1_bin, right1_bin)
  }else if (bin_method == 'M1'){
    bins <- data.frame(bin = left1_bin + right1_bin)
  }else if (bin_method == 'M12'){
    bins <- M12
  }else if (bin_method == 'M24'){
    bins <- M24
  }

  return(bins)

}


#' @title Normalizes, bins and transforms counts
#'
#' @description Normalizes counts by library size,
#' bin using \code{MILLIPEDE} binning method and then take \sQuote{asinh} or \sQuote{log2} transform.
#' @param count_matrix DNase or ATAC-seq read counts matrix.
#' @param idxstats_file The \code{idxstats} file (generated by \code{samtools}).
#' @param ref_size Scale to DNase-seq or ATAC-seq reference library size.
#' (Default: 1e8 for DNase-seq and 5e7 for ATAC-seq).
#' @param bin_method \code{MILLIPEDE} binning scheme (Default: \sQuote{M5}).
#' @param transform Type of transformation for DNase or ATAC counts.
#' Options: \sQuote{asinh}, \sQuote{log2}, \sQuote{sqrt}, \sQuote{none}.
#' @return A data frame of normalized, binned and transformed counts.
#' @export
#' @examples
#' \dontrun{
#' # Normalizes counts by scaling to a library size of 100 million reads,
#' # and bin counts using MILLIPEDE M5 binning method, and then takes
#' # asinh transformation on the binned counts.
#' bins <- normalize_bin_transform_counts(count_matrix,
#'                                        idxstats_file,
#'                                        ref_size = 1e8,
#'                                        bin_method = 'M5',
#'                                        transform = 'asinh')
#' }
normalize_bin_transform_counts <- function(count_matrix,
                                           idxstats_file,
                                           ref_size=1e8,
                                           bin_method = c('M5','M24','M12','M3','M2','M1'),
                                           transform = c('asinh','log2', 'sqrt', 'none')) {

  bin_method <- match.arg(bin_method)
  transform <- match.arg(transform)

  # Normalize (scale) read counts
  normalized_counts <- normalize_counts(count_matrix, idxstats_file, ref_size)

  ## MILLIPEDE binning and transform counts
  bins <- bin_transform_counts(normalized_counts, bin_method, transform)

  return(bins)

}


#' @title Normalizes read counts
#'
#' @description Normalizes DNase or ATAC-seq read counts by library sizes.
#' It first obtains the total mapped reads from the current sample, and
#' then scales the read counts for the current data to a reference library size.
#' @param counts DNase or ATAC-seq read counts matrix
#' @param idxstats_file The \code{idxstats} file generated by \code{samtools}.
#' @param ref_size Normalize to reference library size.
#' (Default: 1e8 for DNase-seq and 5e7 for ATAC-seq).
#' @return A matrix of normalize read counts.
#' @export
normalize_counts <- function(counts,
                             idxstats_file,
                             ref_size=1e8) {

  # Count total mapped reads (chr1:22)
  total_readsMapped <- get_total_reads(idxstats_file, select_chr = TRUE)

  # Normalize (scale) read counts
  cat('Normalize (scale) reads library to', ref_size / 1e6, 'million reads. \n')
  scaling_factor <- ref_size / total_readsMapped
  normalized_counts <- counts * scaling_factor

  return(normalized_counts)

}


#' @title Bins and transforms count matrix
#'
#' @description Bins DNase or ATAC counts
#' using \code{MILLIPEDE} binning and then
#' take \sQuote{sqrt} or \sQuote{log2} transform
#' @param counts DNase or ATAC-seq read counts matrix
#' @param bin_method \code{MILLIPEDE} binning scheme (Default: \sQuote{M5}).
#' @param transform Type of transformation for DNase or ATAC counts.
#' Options: \sQuote{asinh}, \sQuote{log2}, \sQuote{sqrt}, \sQuote{none}.
#' @return A data frame of binned and transformed counts.
#' @export
bin_transform_counts <- function(counts,
                                 bin_method = c('M5','M24','M12','M3','M2','M1'),
                                 transform = c('asinh','log2', 'sqrt', 'none')) {

  bin_method <- match.arg(bin_method)
  transform <- match.arg(transform)

  cat('Perform', bin_method, 'binning... \n')
  bins <- millipede_binning(counts, bin_method)

  if (transform == 'asinh') {
    cat('Perform asinh transform. \n')
    bins <- asinh(bins)
  } else if (transform == 'log2') {
    cat('Perform log2 transform. \n')
    bins <- log2(bins+1)
  } else if (transform == 'sqrt') {
    cat('Perform sqrt transform. \n')
    bins <- sqrt(bins)
  }

  return(bins)

}


#' @title Merges DNase or ATAC-seq counts from multiple replicates,
#' then normalizes merged counts.
#' @description Merges DNase or ATAC-seq read counts from multiple replicate samples,
#' and normalizes the read counts by scaling to the reference library size.
#'
#' @param counts_files DNase or ATAC-seq read counts matrix files.
#' @param idxstats_files The \code{idxstats} files generated by \code{samtools}.
#' @param ref_size Scales to reference library size.
#' (Default: 1e8 for DNase-seq and 5e7 for ATAC-seq).
#' @return A matrix of merged and normalized counts.
#' @export
merge_normalize_counts <- function(counts_files, idxstats_files, ref_size = 1e8){

  ## Loads raw counts and merges the replicates
  for (i in 1:length(counts_files)) {
    count_matrix <- readRDS(counts_files[i])
    if ( i == 1 ) {
      total_countmatrix <- count_matrix
    }else{
      total_countmatrix <- total_countmatrix + count_matrix
    }
  }

  # Counts total mapped reads
  total_readsMapped <- sum(sapply(idxstats_files, get_total_reads, select_chr = TRUE))

  # Normalizes (scales) read counts
  cat('Normalize (scale) to', ref_size / 1e6, 'million reads. \n')
  scaling_factor <- ref_size / total_readsMapped
  normalized_countmatrix <- total_countmatrix * scaling_factor

  return(normalized_countmatrix)

}

#' @title Merges DNase or ATAC-seq counts from multiple replicates,
#' then normalizes, bin and transform the merged counts
#' @description Merges DNase-seq or ATAC-seq read counts from
#' multiple replicate samples,
#' normalizes the read counts by scaling to the reference library size,
#' then bins and transforms the merged counts.
#'
#' @param counts_files DNase or ATAC-seq read counts matrix files.
#' @param idxstats_files The \code{idxstats} files generated by samtools.
#' @param ref_size Scale to DNase-seq or ATAC-seq reference library size.
#' (Default: 1e8 for DNase-seq and 5e7 for ATAC-seq).
#' @param bin_method \code{MILLIPEDE} binning scheme (default: \sQuote{M5}).
#' @param transform Type of transformation for DNase or ATAC counts.
#' Options: \sQuote{asinh}, \sQuote{log2}, \sQuote{sqrt}, \sQuote{none}.
#' @return A data frame of merged, normalized, binned, and transformed counts.
#' @export
merge_normalize_bin_transform_counts <- function(counts_files,
                                                 idxstats_files,
                                                 ref_size = 1e8,
                                                 bin_method = c('M5','M24','M12','M3','M2','M1'),
                                                 transform = c('asinh','log2', 'sqrt', 'none')){


  # Merge counts from multiple replicates, then normalize merged counts
  normalized_countmatrix <- merge_normalize_counts(counts_files, idxstats_files, ref_size)

  # MILLIPEDE binning and transform counts
  bins <- bin_transform_counts(normalized_countmatrix, bin_method, transform)

  return(bins)

}


# Converts bedGraph to BigWig format
bedGraphToBigWig <- function(bedgraph_file,
                             chrom_size_file,
                             bedGraphToBigWig_path='bedGraphToBigWig'){

  # Checking input arguments
  if ( Sys.which(bedGraphToBigWig_path) == '' )
    stop( 'bedGraphToBigWig could not be executed. Please install bedGraphToBigWig and set bedGraphToBigWig_path.' )

  # Convert bedGraph to BigWig format
  outdir <- dirname(bedgraph_file)
  outname <- tools::file_path_sans_ext(basename(bedgraph_file))
  bw_file <- file.path(outdir, paste0(outname, '.bw'))

  cmd <- paste(bedGraphToBigWig_path, bedgraph_file, chrom_size_file, bw_file)
  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

}

