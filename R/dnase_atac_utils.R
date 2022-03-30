
#' @title Count DNase or ATAC cuts around candidate sites
#' @description Count DNase or ATAC cuts around candidate sites on both strands.
#' Then combine the counts into one matrix, with the first half of the columns
#' representing the read counts on the forward strand,
#' and the second half of the columns representing the read counts
#' on the reverse strand with respect to the motif match location and strand.
#' @param sites A data frame containing the candidate sites.
#' @param bam_file BAM file of DNase or ATAC-seq reads alignment
#' @param shift_ATAC Logical. When \code{shift_ATAC = TRUE},
#' it shifts reads aligned to the + strand by +4 bp,
#' and shifts reads aligned to the - strand by -5 bp.
#' @param site_0base Logical. If TRUE,
#' the coordinates of sites are 0-based as in BED format. Otherwise,
#' the coordinates of sites are 1-based.
#' @param ncores Number of cores available to parallel computing.
#' @return A count matrix. The first half of the columns
#' are the read counts on the forward strand, and the second half of the
#' columns are the read counts on the reverse strand.
#' @export
#' @examples
#' # Count DNase or ATAC cuts around candidate sites
#' sites_counts.mat <- count_sites_cuts(sites,
#'                                      bam_file = "ATAC_example.bam",
#'                                      shift_ATAC = TRUE,
#'                                      site_0base = TRUE,
#'                                      ncores = 10)
count_sites_cuts <- function(sites,
                             bam_file,
                             shift_ATAC = FALSE,
                             site_0base = TRUE,
                             ncores = 10){

  # read cuts from BAM file
  cuts.gr <- read_bam_cuts(bam_file, shift_ATAC, return_type = 'cuts')

  if(site_0base)
    sites$start <- sites$start  + 1

  sites.gr <- makeGRangesFromDataFrame(sites)

  cuts.gr <- subsetByOverlaps(cuts.gr, sites.gr, ignore.strand = TRUE)

  n_available_cores <- detectCores(logical = FALSE) - 1
  if(missing(ncores)){
    ncores <- n_available_cores
  }else{
    ncores <- min(n_available_cores, ncores)
  }

  cat(sprintf('Counting cuts using %d cores...\n', ncores))
  registerDoParallel(cores=ncores)

  fwd_counts.mat <- foreach(i=1:length(sites.gr), .combine=rbind) %dopar% {
    countOverlaps(query = GPos(sites.gr[i]), subject = cuts.gr[strand(cuts.gr) == '+'],
                  ignore.strand = TRUE)
  }

  rev_counts.mat <- foreach(i=1:length(sites.gr), .combine=rbind) %dopar% {
    countOverlaps(query = GPos(sites.gr[i]), subject = cuts.gr[strand(cuts.gr) == '-'],
                  ignore.strand = TRUE)
  }

  sites_counts.l <- list(fwd = fwd_counts.mat, rev = rev_counts.mat)
  neg_strand <- (sites$strand == '-')
  sites_counts.l$fwd[neg_strand, ] <- t(apply(rev_counts.mat[neg_strand, ], 1, rev))
  sites_counts.l$rev[neg_strand, ] <- t(apply(fwd_counts.mat[neg_strand, ], 1, rev))

  sites_counts.mat <- as.matrix(cbind(sites_counts.l$fwd, sites_counts.l$rev))
  rownames(sites_counts.mat) <- sites$name
  colnames(sites_counts.mat) <- c(paste0('F', c(1:ncol(fwd_counts.mat))),
                                  paste0('R', c(1:ncol(rev_counts.mat))))

  return(sites_counts.mat)
}


#' @title Get count matrices around candidate binding sites
#' @description Extract counts around candidate binding sites on both strands
#' from the genome counts (BigWig files prepared using \code{count_genome_cuts}).
#' It utilizes the \code{extract bed} function from the \code{bwtool} software
#' to extract the read counts.
#' Then combine the counts into one matrix, with the first half of the columns
#' representing the read counts on the forward strand,
#' and the second half of the columns representing the read counts
#' on the reverse strand.
#' @param sites A data frame containing the candidate sites.
#' @param genomecount_dir Directory for genome counts.
#' @param genomecount_name File prefix for genome counts.
#' @param tmpdir Temporary directory to save intermediate files.
#' @param bwtool_path Path to \code{bwtool} executable.
#' @return A count matrix. The first half of the columns
#' are the read counts on the forward strand, and the second half of the
#' columns are the read counts on the reverse strand.
#' @importFrom data.table fwrite
#' @export
#' @examples
#' # Get ATAC-seq count matrices around candidate sites
#' sites_counts.mat <- get_sites_counts_bwtool(sites,
#'                                             bam_file='K562.ATAC.bam',
#'                                             chrom_size_file,
#'                                             shift_ATAC=TRUE,
#'                                             genomecount_dir='processed_data',
#'                                             genomecount_name='K562.ATAC',
#'                                             bedGraphToBigWig_path='bedGraphToBigWig',
#'                                             bwtool_path='bwtool')
get_sites_counts_bwtool <- function(sites,
                                    bam_file,
                                    chrom_size_file,
                                    shift_ATAC=FALSE,
                                    genomecount_dir,
                                    genomecount_name,
                                    tmpdir=genomecount_dir,
                                    bedGraphToBigWig_path='bedGraphToBigWig',
                                    bwtool_path='bwtool') {

  if ( Sys.which(bwtool_path) == '' ) {
    stop( 'bwtool could not be executed. Please install bwtool and set bwtool_path.' )
  }

  genome_fwd_count_file <- file.path(genomecount_dir, paste0(genomecount_name, '.fwd.genomecounts.bw'))
  genome_rev_count_file <- file.path(genomecount_dir, paste0(genomecount_name, '.rev.genomecounts.bw'))

  if( !(file.exists(genome_fwd_count_file) & file.exists(genome_rev_count_file)) ){
    count_genome_cuts_bw(bam_file, chrom_size_file, shift_ATAC,
                         genomecount_dir, genomecount_name)
  }

  fwd_matrix_file <- tempfile(pattern = 'sitescounts.', tmpdir = tmpdir, fileext = '.fwd.matrix')
  rev_matrix_file <- tempfile(pattern = 'sitescounts.', tmpdir = tmpdir, fileext = '.rev.matrix')

  sites_file <- tempfile(pattern = 'sites.', tmpdir = tmpdir, fileext = '.txt')
  fwrite(sites[,1:4], sites_file, sep = '\t', col.names = FALSE, scipen = 999)

  cat('Extract counts around candidate sites ... \n')
  cmd <- paste(bwtool_path, 'extract bed', sites_file,
               genome_fwd_count_file, fwd_matrix_file, '-fill=0 -decimals=0 -tabs')
  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

  cmd <- paste(bwtool_path, 'extract bed', sites_file,
               genome_rev_count_file, rev_matrix_file, '-fill=0 -decimals=0 -tabs')
  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

  # Flip the counts generated from bwtool for motifs on the reverse strand
  sites_counts.l <- flip_rev_strand_counts(sites, fwd_matrix_file, rev_matrix_file)

  # Combine counts on both strands
  sites_counts.mat <- as.matrix(cbind(sites_counts.l$fwd, sites_counts.l$rev))
  rownames(sites_counts.mat) <- sites$name

  unlink(c(sites_file, fwd_matrix_file, rev_matrix_file))
  return(sites_counts.mat)
}

#' @title Count DNase-seq or ATAC-seq cuts along the genome
#' @description Count genomic cleavage (5' end) from DNase-seq or
#' ATAC-seq BAM alignment files.
#' For ATAC-seq, when \code{shift_ATAC = TRUE}, it shifts reads aligned to
#' the + strand by +4 bp, and shifts reads aligned to the - strand
#' by -5 bp (Buenrostro et al. 2013).
#' This function treats paired end reads as independent.
#' If you have paired-end ATAC-seq reads and want to select the fragments
#' within NFR regions, please refer to the
#' \code{count_ATAC_genome_cuts_readpairs()} function.
#' @param bam_file Sorted BAM file.
#' @param chrom_size_file File of genome sizes by chromosomes.
#' @param shift_ATAC Logical. When \code{shift_ATAC = TRUE},
#' it shifts reads aligned to the + strand by +4 bp,
#' and shifts reads aligned to the - strand by -5 bp.
#' @param outdir Output directory (default: use the directory of \code{bam_file}).
#' @param outname Output prefix (default: use the prefix of \code{bam_file}).
#' @param bedGraphToBigWig_path Path to UCSC \code{bedGraphToBigWig} executable.#'
#' @importFrom GenomicAlignments readGAlignments
#' @import GenomicRanges
#' @importFrom data.table fwrite
#' @export
count_genome_cuts_bw <- function(bam_file,
                                 chrom_size_file,
                                 shift_ATAC = FALSE,
                                 outdir = dirname(bam_file),
                                 outname,
                                 bedGraphToBigWig_path='bedGraphToBigWig'){

  cat('Counting genome cuts for', bam_file, '...\n')
  coverage.gr <- read_bam_cuts(bam_file, shift_ATAC, return_type = 'coverage')
  pos_coverage.gr <- coverage.gr$pos
  neg_coverage.gr <- coverage.gr$neg

  if(!dir.exists(outdir))
    dir.create(outdir)

  if(missing(outname))
    outname <- tools::file_path_sans_ext(basename(bam_file))

  # + strand counts
  cat('Save genome counts for + strand ...\n')
  pos_coverage.df <- data.frame(chr = seqnames(pos_coverage.gr),
                                start = start(pos_coverage.gr) - 1,
                                end = end(pos_coverage.gr),
                                count = pos_coverage.gr$score)
  pos_coverage.df <- pos_coverage.df[pos_coverage.df$count > 0, ]
  bedgraph_file <- file.path(outdir, paste0(outname, '.fwd.genomecounts.bedGraph'))
  data.table::fwrite(pos_coverage.df, bedgraph_file,
                     sep='\t', col.names = FALSE, scipen = 999)
  bedGraphToBigWig(bedgraph_file, chrom_size_file, bedGraphToBigWig_path)
  unlink(bedgraph_file)

  # - strand counts
  cat('Save genome counts for - strand ...\n')
  neg_coverage.df <- data.frame(chr = seqnames(neg_coverage.gr),
                                start = start(neg_coverage.gr) - 1,
                                end = end(neg_coverage.gr),
                                count = neg_coverage.gr$score)
  neg_coverage.df <- neg_coverage.df[neg_coverage.df$count > 0, ]
  bedgraph_file <- file.path(outdir, paste0(outname, '.rev.genomecounts.bedGraph'))
  data.table::fwrite(neg_coverage.df, bedgraph_file,
                     sep='\t', col.names = FALSE, scipen = 999)
  bedGraphToBigWig(bedgraph_file, chrom_size_file, bedGraphToBigWig_path)
  unlink(bedgraph_file)

}

#' @title Count DNase-seq or ATAC-seq cuts along the genome
#' @description Count genomic cleavage (5' end) from DNase-seq or
#' ATAC-seq BAM alignment files using \code{bedtools}
#' For ATAC-seq, when \code{shift_ATAC = TRUE}, it shifts reads aligned to
#' the + strand by +4 bp, and shifts reads aligned to the - strand
#' by -5 bp (Buenrostro et al. 2013).
#' This function treats paired end reads as independent.
#' If you have paired-end ATAC-seq reads and want to select the fragments
#' within NFR regions, please refer to the
#' \code{count_ATAC_genome_cuts_readpairs()} function.
#' @param bam_file Sorted BAM file.
#' @param chrom_size_file File of genome sizes by chromosomes.
#' @param shift_ATAC Logical. When \code{shift_ATAC = TRUE},
#' it shifts reads aligned to the + strand by +4 bp,
#' and shifts reads aligned to the - strand by -5 bp.
#' @param outdir Output directory (default: use the directory of \code{bam_file}).
#' @param outname Output prefix (default: use the prefix of \code{bam_file}).
#' @param bedtools_path Path to \code{bedtools} executable.
#' @param bedGraphToBigWig_path Path to UCSC \code{bedGraphToBigWig} executable.
#' @importFrom GenomicAlignments readGAlignments
#' @import GenomicRanges
#' @importFrom data.table fwrite
#' @export
count_genome_cuts_bedtools <- function(bam_file,
                                       chrom_size_file,
                                       shift_ATAC=FALSE,
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

  if(missing(outdir))
    outdir <- dirname(bam_file)

  if(!dir.exists(outdir))
    dir.create(outdir)

  if(missing(outname))
    outname <- tools::file_path_sans_ext(basename(bam_file))

  for (strand in c('+', '-')) {
    if (strand == '+') {
      bedgraph_file <- file.path(outdir, paste0(outname, '.fwd.genomecounts.bedGraph'))
    } else {
      bedgraph_file <- file.path(outdir, paste0(outname, '.rev.genomecounts.bedGraph'))
    }

    # Count 5' end cleavage and save as a .bedGraph file
    cat('Counting genome cleavage for', bam_file, 'on', strand, 'strand...\n')
    cmd <- paste(bedtools_path, 'genomecov -bg -5', '-strand', strand,
                 '-ibam', bam_file, '>', bedgraph_file)
    if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

    cat('Sorting counted genome cleavage on', strand, 'strand...\n')
    genome_counts <- data.table::fread(bedgraph_file)
    genome_counts <- genome_counts[with(genome_counts, order(V1, V2)), ]
    # Shift ATAC-seq reads to get the centers of Tn5 binding positions
    # Shift +4 bp for + strand and -5 bp for - strand.
    if (shift_ATAC) {
      cat('Shifting ATAC-seq reads ...\n')
      if (strand == '+') {
        genome_counts[,c(2:3)] <- genome_counts[,c(2:3)] + 4
      } else {
        genome_counts[,c(2:3)] <- genome_counts[,c(2:3)] - 5
      }
    }

    data.table::fwrite(genome_counts, bedgraph_file,
                       sep='\t', col.names = FALSE, scipen = 999)

    # save to bigwig file
    bedGraphToBigWig(bedgraph_file, chrom_size_file, bedGraphToBigWig_path)
    unlink(bedgraph_file)
  }

}



#' @title Flip the counts generated from \code{bwtool} for motifs
#' on the reverse strand
#'
#' @param sites A data frame containing candidate sites
#' @param fwd_matrix_file File for count matrix on the forward strand
#' @param rev_matrix_file File for count matrix on the reverse strand
#' @importFrom  data.table fread fwrite
flip_rev_strand_counts <- function(sites, fwd_matrix_file, rev_matrix_file) {

  fwd_count <- as.data.frame(fread(fwd_matrix_file))
  rev_count <- as.data.frame(fread(rev_matrix_file))

  if(sum(sites[,2] != fwd_count[,2]) != 0) {
    stop('Sites do not match!')
  }

  # Extract the count values
  fwd_count.m <- as.matrix(fwd_count[, -c(1:5)])
  colnames(fwd_count.m) <- paste0('F', c(1:ncol(fwd_count.m)))
  rev_count.m <- as.matrix(rev_count[, -c(1:5)])
  colnames(rev_count.m) <- paste0('R', c(1:ncol(rev_count.m)))

  # For motifs match to the - strand, flip the fwd and rev counts, and reverse the counts
  sites_counts.l <- list(fwd = fwd_count.m, rev = rev_count.m)
  minusStrand <- (sites$strand == '-')
  sites_counts.l$fwd[minusStrand, ] <- t(apply(rev_count.m[minusStrand, ], 1, rev))
  sites_counts.l$rev[minusStrand, ] <- t(apply(fwd_count.m[minusStrand, ], 1, rev))

  fwd_count <- cbind(sites[,c(1:3,6)], sites_counts.l$fwd)
  rev_count <- cbind(sites[,c(1:3,6)], sites_counts.l$rev)

  fwrite(fwd_count, fwd_matrix_file, sep = ' ', scipen = 999)
  fwrite(rev_count, rev_matrix_file, sep = ' ', scipen = 999)

  return(sites_counts.l)
}


#' Read DNase or ATAC cuts or coverage of cuts from BAM file
#'
#' @param bam_file BAM file of DNase or ATAC-seq reads alignment
#' @param shift_ATAC Logical. When \code{shift_ATAC = TRUE},
#' it shifts reads aligned to the + strand by +4 bp,
#' and shifts reads aligned to the - strand by -5 bp.
#' @param return_type Options for the returned data:
#' \dQuote{cuts} or \dQuote{coverage}.
#' @return A GRange object of cuts or coverage
read_bam_cuts <- function(bam_file,
                          shift_ATAC = FALSE,
                          return_type = c('cuts', 'coverage')){

  return_type <- match.arg(return_type)
  cat('Reading genome cuts for', bam_file, '...\n')

  # Read ATAC-seq alignments
  reads <- GenomicAlignments::readGAlignments(bam_file)

  # Extract 5' end position and shift based on strand
  if(shift_ATAC){
    cat('Shift ATAC-seq reads... \n')
    cuts <- GenomicRanges::resize(GenomicRanges::granges(reads), fix = 'start', 1)
    pos_cuts <- GenomicRanges::shift(GenomicRanges::granges(cuts[strand(cuts) == '+']), 4)
    neg_cuts <- GenomicRanges::shift(GenomicRanges::granges(cuts[strand(cuts) == '-']), -5)
    cuts <- c(pos_cuts, neg_cuts)
  }else{
    cuts <- GenomicRanges::resize(GenomicRanges::granges(reads), fix = 'start', 1)
  }

  if(return_type == 'coverage'){
    pos_coverage.gr <- GRanges(coverage(cuts[strand(cuts) == '+']), strand = '+')
    pos_coverage.gr <- sort(pos_coverage.gr)
    neg_coverage.gr <- GRanges(coverage(cuts[strand(cuts) == '-']), strand = '-')
    neg_coverage.gr <- sort(neg_coverage.gr)

    cuts_coverage.gr <- GRangesList(pos = pos_coverage.gr,
                                    neg = neg_coverage.gr)
    return(cuts_coverage.gr)
  }else if(return_type == 'cuts'){
    return(cuts)
  }

}

#' @title Perform MILLIPEDE binning on count matrix
#' @description Perform binning using different MILLIPEDE binning schemes
#' (M5, M24, M12, M3, M2, M1) on the input count matrix.
#' @param counts DNase-seq or ATAC-seq read counts matrix,
#' rows are candidate sites,
#' columns are DNase or ATAC counts with 100bp flanks around motifs
#' on the forward and reverse strands.
#' @param bin_method MILLIPEDE binning scheme. Options:
#' \dQuote{M5} (default), \dQuote{M24}, \dQuote{M12}, \dQuote{M3}, \dQuote{M2},
#' and \dQuote{M1}.
#' @param combine_strands Method to combine counts on both strands from M24 bins
#' to M12 bins: vertical' (default) or 'motif'.
#' @return A list containing binning results (data frames) using different
#' binning schemes.
#' @export
#' @examples
#'
#' # Perform MILLIPEDE binning with different binning schemes.
#'
#' # M5 binning
#'  M5_bins <- millipede_binning(counts, bin_method = 'M5')
#'
#' # M12 binning
#'  M12_bins <- millipede_binning(counts, bin_method = 'M12')
#'
#' # M3 binning
#'  M3_bins <- millipede_binning(counts, bin_method = 'M3')
#'
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


#' @title Normalize, bin and transform counts
#'
#' @description Normalize counts by library size,
#' bin using MILLIPEDE binning method and then take asinh or log2 transform
#'
#' @param count_matrix DNase or ATAC-seq read counts matrix.
#' @param idxstats_file The \code{idxstats} file (generated by \code{samtools}).
#' @param ref_size Scale to DNase-seq or ATAC-seq reference library size (Default: 1e8).
#' @param bin_method MILLIPEDE binning scheme (Default: \dQuote{M5}).
#' @param transform Type of transformation for DNase or ATAC counts.
#' Options: \dQuote{asinh}, \dQuote{log2}, \dQuote{sqrt}, \dQuote{none}.
#' @return A data frame of normalized, binned and transformed counts.
#' @export
#' @examples
#'
#' # Normalize counts by scaling to a library size of 100 million reads,
#' # and bin counts using MILLIPEDE M5 binning method, and then take
#' # asinh transformation on the binned counts.
#' bins <- normalize_bin_transform_counts(count_matrix,
#'                                        idxstats_file,
#'                                        ref_size = 1e8,
#'                                        bin_method = 'M5',
#'                                        transform = 'asinh')
#'
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


#' @title Normalize read counts
#'
#' @description Normalize DNase or ATAC-seq read counts by library sizes.
#' It first obtain the total mapped reads from the current sample, and
#' then scales the read counts for the current data to a reference library size.
#' @param counts DNase or ATAC-seq read counts matrix
#' @param idxstats_file The \code{idxstats} file generated by samtools.
#' @param ref_size Normalize to reference library size (default: 1e8).
#' @return A matrix of normalize read counts.
#'
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


#' @title Bin and transform count matrix
#'
#' @description Binning DNase or ATAC count matrix
#' using MILLIPEDE binning and then take asinh or log2 transform
#' @param counts DNase or ATAC-seq read counts matrix
#' @param bin_method MILLIPEDE binning scheme (Default: \dQuote{M5}).
#' @param transform Type of transformation for DNase or ATAC counts.
#' Options: \dQuote{asinh}, \dQuote{log2}, \dQuote{sqrt}, \dQuote{none}.
#' @return A data frame of binned and transformed counts.
#'
bin_transform_counts <- function(counts,
                                 bin_method = c('M5','M24','M12','M3','M2','M1'),
                                 transform = c('asinh','log2', 'sqrt', 'none')) {

  bin_method <- match.arg(bin_method)
  transform <- match.arg(transform)

  cat('Perform MILLIPEDE', bin_method, 'binning... \n')
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


#' @title Merge DNase or ATAC-seq counts from multiple replicates,
#' then normalize merged counts.
#' @description Merge DNase or ATAC-seq read counts from multiple replicate samples,
#' and normalize the read counts by scaling to the reference library size.
#'
#' @param counts_files DNase or ATAC-seq read counts matrix files.
#' @param idxstats_files The idxstats files generated by \code{samtools}.
#' @param ref_size Scale to reference library size (default: 1e8).
#' @return A matrix of merged and normalized counts.
#' @export
#'
merge_normalize_counts <- function(counts_files, idxstats_files, ref_size = 1e8){

  ## Load raw counts and merge the replicates
  for (i in 1:length(counts_files)) {
    count_matrix <- readRDS(counts_files[i])
    if ( i == 1 ) {
      total_countmatrix <- count_matrix
    }else{
      total_countmatrix <- total_countmatrix + count_matrix
    }
  }

  # Count total mapped reads
  total_readsMapped <- sum(sapply(idxstats_files, get_total_reads, select_chr = TRUE))

  # Normalize (scale) read counts
  cat('Normalize (scale) to', ref_size / 1e6, 'million reads. \n')
  scaling_factor <- ref_size / total_readsMapped
  normalized_countmatrix <- total_countmatrix * scaling_factor

  return(normalized_countmatrix)

}

#' @title Merge DNase or ATAC-seq counts from multiple replicates,
#' then normalize, bin and transform the merged counts
#' @description Merge DNase-seq or ATAC-seq read counts from
#' multiple replicate samples,
#' normalize the read counts by scaling to the reference library size,
#' then bin and transform the merged counts.
#'
#' @param counts_files DNase or ATAC-seq read counts matrix files.
#' @param idxstats_files The idxstats files generated by samtools.
#' @param ref_size Scale to DNase-seq or ATAC-seq reference library size
#' (default: 1e8).
#' @param bin_method MILLIPEDE binning scheme (default: \dQuote{M5}).
#' @param transform Type of transformation for DNase or ATAC counts.
#' Options: \dQuote{asinh}, \dQuote{log2}, \dQuote{sqrt}, \dQuote{none}.
#' @return A data frame of merged, normalized, binned, and transformed counts.
#' @export
#'
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


#' @title Count genome coverage of shifted ATAC 5' cuts for ATAC-seq paired-end reads
#' @param bam_file BAM file of ATAC-seq paired-end reads
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @param outdir Output directory. Default: the directory of the \code{bam_file}
#' @param plot_fragsize Logical, if TRUE, plot the distribution of fragment sizes
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom rtracklayer export
#' @import GenomicRanges
#' @importFrom data.table fwrite
count_ATAC_genome_cuts_readpairs <- function(bam_file,
                                               select_NFR_fragments = FALSE,
                                               shift_ATAC = FALSE,
                                               write_bedgraph = FALSE,
                                               outdir = dirname(bam_file),
                                               outname){

  cat('Counting genome cuts for', bam_file, '...\n')

  # Load only properly paired reads
  flags <- Rsamtools::scanBamFlag(isProperPair = TRUE)
  param <- Rsamtools::ScanBamParam(flag = flags,
                                   what = c('qname', 'flag', 'mapq', 'isize'))
  # Read ATAC-seq alignment pairs
  readpairs <- GenomicAlignments::readGAlignmentPairs(bam_file, param = param)

  # Separate read pairs
  read1 <- first(readpairs)
  read2 <- second(readpairs)

  if(select_NFR_fragments){
    # Fragment insert sizes
    insert_sizes <- abs(elementMetadata(read1)$isize)
    # Select fragments < 100bp (within nucleosome free regions)
    cat('Select ATAC-seq fragments < 100bp (within nucleosome free regions) \n')
    selected_readpairs <- readpairs[insert_sizes < 100, ]
    selected_bam_file <- file.path(outdir, gsub('\\.bam', '_NFR.bam', basename(bam_file)))

    cat('Save filtered bam file to:', selected_bam_file, '\n')
    rtracklayer::export(selected_readpairs, selected_bam_file, format = 'bam')

  }

  # Extract 5' end position and shift based on strand
  if(shift_ATAC){
    cat('Shift ATAC-seq reads... \n')
    read1_5ends <- resize(granges(read1), fix = 'start', 1)
    read1_pos_cuts <- shift(granges(read1_5ends[strand(read1_5ends) == '+']), 4)
    read1_neg_cuts <- shift(granges(read1_5ends[strand(read1_5ends) == '-']), -5)

    read2_5ends <- resize(granges(read2), fix = 'start', 1)
    read2_pos_cuts <- shift(granges(read2_5ends[strand(read2_5ends) == '+']), 4)
    read2_neg_cuts <- shift(granges(read2_5ends[strand(read2_5ends) == '-']), -5)

    combined_cuts <- c(read1_pos_cuts, read1_neg_cuts,
                       read2_pos_cuts, read2_neg_cuts)
  }else{
    read1_5ends <- resize(granges(read1), fix = 'start', 1)
    read2_5ends <- resize(granges(read2), fix = 'start', 1)
    combined_cuts <- c(read1_5ends, read2_5ends)
  }

  pos_coverage.gr <- GRanges(coverage(combined_cuts[strand(combined_cuts) == '+']), strand = '+')
  pos_coverage.gr <- sort(pos_coverage.gr)
  neg_coverage.gr <- GRanges(coverage(combined_cuts[strand(combined_cuts) == '-']), strand = '-')
  neg_coverage.gr <- sort(neg_coverage.gr)
  cuts_coverage.gr <- GRangesList(pos = pos_coverage.gr,
                                  neg = neg_coverage.gr)

  if(write_bedgraph){

    if(!dir.exists(outdir))
      dir.create(outdir)

    if(missing(outname))
      outname <- tools::file_path_sans_ext(basename(bam_file))

    if(shift_ATAC)
      outname <- paste0(outname, '.shifted')
    # + strand counts
    pos_coverage.df <- data.frame(chr = seqnames(pos_coverage.gr),
                                  start = start(pos_coverage.gr) - 1,
                                  end = end(pos_coverage.gr),
                                  count = pos_coverage.gr$score)
    pos_coverage.df <- pos_coverage.df[pos_coverage.df$count > 0, ]
    bedgraph_file <- file.path(outdir, paste0(outname, '.pos.genomecounts.bedGraph'))
    data.table::fwrite(pos_coverage.df, bedgraph_file,
                       sep='\t', col.names = FALSE, scipen = 999)
    # - strand counts
    neg_coverage.df <- data.frame(chr = seqnames(neg_coverage.gr),
                                  start = start(neg_coverage.gr) - 1,
                                  end = end(neg_coverage.gr),
                                  count = neg_coverage.gr$score)
    neg_coverage.df <- neg_coverage.df[neg_coverage.df$count > 0, ]
    bedgraph_file <- file.path(outdir, paste0(outname, '.neg.genomecounts.bedGraph'))
    data.table::fwrite(neg_coverage.df, bedgraph_file,
                       sep='\t', col.names = FALSE, scipen = 999)
  }

  return(cuts_coverage.gr)
}

# Convert bedGraph to BigWig format
bedGraphToBigWig <- function(bedgraph_file,
                             chrom_size_file,
                             bedGraphToBigWig_path){

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



