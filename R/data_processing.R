
#' @title Sort and index the BAM file, and then retrieve and print stats in the index file.
#'
#' @param bam_file Input BAM file.
#' @param dir_output Output directory.
#' @param sort logical. If TRUE, sort (and index) the BAM file.
#' @param stats logical. If TRUE, retrieve and print stats in the index file.
#' @param path_samtools Path to samtools executable.
#'
#' @export
bam_sort_index_stats <- function(bam_file, dir_output=NA, sort=TRUE, stats=TRUE,
                                 path_samtools='samtools') {

  if ( system( paste(path_samtools,'--help') , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( 'ERROR: samtools could not be executed, set path_samtools\n' , sep='', file=stderr() )
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

  if ( system( paste(path_samtools,'--help') , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( 'ERROR: samtools could not be executed, set path_samtools\n' , sep='', file=stderr() )
    cleanup()
    q()
  }

  if ( system( paste(path_bedtools,'--help') , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( 'ERROR: bedtools could not be executed, set path_bedtools\n' , sep='', file=stderr() )
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

    if ( system( paste(path_bedGraphToBigWig,'--help') , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
      cat( 'ERROR: bedGraphToBigWig could not be executed, set path_bedGraphToBigWig\n' , sep='', file=stderr() )
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

  cat('Flipping tagcounts for motifs on the reverse strand ... \n')

  fwd_count.df <- read.table(file_dnase_matrix_fwd)
  rev_count.df <- read.table(file_dnase_matrix_rev)
  sites.df <- read.table(file_sites)

  if(sum(sites.df[,2] != fwd_count.df[,2]) != 0) {
    stop('Sites do not match!')

  } else {

    # Extract the count values
    fwd_count.m <- fwd_count.df[, 6:ncol(fwd_count.df)]
    rev_count.m <- rev_count.df[, 6:ncol(rev_count.df)]

    fwd_output <- fwd_count.m
    rev_output <- rev_count.m

    # For motifs match to the minus strand, flip the fwd and rev tagcounts, and reverse the counts
    idx_minusStrand <- which(sites.df[,6] == '-')

    fwd_output[idx_minusStrand, ] <- t(apply(rev_count.m[idx_minusStrand, ], 1, rev))
    rev_output[idx_minusStrand, ] <- t(apply(fwd_count.m[idx_minusStrand, ], 1, rev))

    # Add the sites info to the first few columns
    fwd_count.df <- cbind(sites.df[,c(1:3,6)], fwd_output)
    rev_count.df <- cbind(sites.df[,c(1:3,6)], rev_output)

    write.table(fwd_count.df, gzfile(file_dnase_matrix_fwd), sep = ' ', quote = F, row.names = F, col.names = F)
    write.table(rev_count.df, gzfile(file_dnase_matrix_rev), sep = ' ', quote = F, row.names = F, col.names = F)

  }
}

#' @title Combine DNase tagcount matrices in both forward and reverse strands
#'
#' @param file_dnase_matrix_fwd filename for DNase tagcount matrix in forward strand
#' @param file_dnase_matrix_rev filename for DNase tagcount matrix in reverse strand
combine_dnase_counts <- function(file_dnase_matrix_fwd, file_dnase_matrix_rev) {
  cat('Combine DNase count matrices in both strands\n')

  DNase_fwd <- read.table(gzfile(file_dnase_matrix_fwd))
  DNase_rev <- read.table(gzfile(file_dnase_matrix_rev))

  DNase_fwd <- DNase_fwd[,5:ncol(DNase_fwd)]
  DNase_rev <- DNase_rev[,5:ncol(DNase_rev)]

  DNase_combined <- cbind(DNase_fwd, DNase_rev)
  return(DNase_combined)
}


#' @title Get the number of total mapped reads from the idxstats file generated by samtools
#'
#' @param file_idxstats idxstats file generated by samtools.
#' @param select_chr whether to select chromosomes or use all chromosomes in the idxstats file.
#' @param chr_list Chromosomes to be included.
#'
#' @return
readstats_total <- function(file_idxstats, select_chr = TRUE, chr_list = paste0('chr', c(1:22))) {
  count.df <- read.table(file_idxstats, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  names(count.df) = c('chr', 'length', 'reads_mapped', 'reads_unmapped')

  if( select_chr == T ){
    total_count <- sum(count.df[which(count.df$chr %in% chr_list), 'reads_mapped'])
  } else {
    total_count <- sum(count.df[, 'reads_mapped'])
  }

  return(total_count)
}

#' @title Count DNase cuts with different MILLIPEDE binning settings.
#'
#' @param DNaseData DNase data matrix, rows are candidate sites,
#' columns are DNase cuts with 100bp flanks around motifs in forward and reverse strands
#' @param combine_strands Method to combine DNase cuts from M24 bins in both strands to M12 bins:
#' 'vertical' (default) or 'motif'.
#'
#' @return
#' @export
sum_cuts_strand <- function(DNaseData, combine_strands='vertical') {

  flank_LF1 = 81:100
  flank_LF2 = flank_LF1 - 20
  flank_LF3 = flank_LF2 - 20
  flank_LF4 = flank_LF3 - 20
  flank_LF5 = flank_LF4 - 20
  flank_RF1 = (ncol(DNaseData)/2-99): (ncol(DNaseData)/2-80)
  flank_RF2 = flank_RF1 + 20
  flank_RF3 = flank_RF2 + 20
  flank_RF4 = flank_RF3 + 20
  flank_RF5 = flank_RF4 + 20

  flank_LR1 = ncol(DNaseData)/2 + flank_LF1
  flank_LR2 = ncol(DNaseData)/2 + flank_LF2
  flank_LR3 = ncol(DNaseData)/2 + flank_LF3
  flank_LR4 = ncol(DNaseData)/2 + flank_LF4
  flank_LR5 = ncol(DNaseData)/2 + flank_LF5
  flank_RR1 = ncol(DNaseData)/2 + flank_RF1
  flank_RR2 = ncol(DNaseData)/2 + flank_RF2
  flank_RR3 = ncol(DNaseData)/2 + flank_RF3
  flank_RR4 = ncol(DNaseData)/2 + flank_RF4
  flank_RR5 = ncol(DNaseData)/2 + flank_RF5

  motif_LF = 101: round(ncol(DNaseData)/4)
  motif_RF = (round(ncol(DNaseData)/4)+1) : (ncol(DNaseData)/2-100)
  motif_LR = ncol(DNaseData)/2 + motif_LF
  motif_RR = ncol(DNaseData)/2 + motif_RF

  M24 = data.frame(
    flank_LF5 = rowSums(DNaseData[,flank_LF5]),
    flank_LF4 = rowSums(DNaseData[,flank_LF4]),
    flank_LF3 = rowSums(DNaseData[,flank_LF3]),
    flank_LF2 = rowSums(DNaseData[,flank_LF2]),
    flank_LF1 = rowSums(DNaseData[,flank_LF1]),

    motif_LF = rowSums(DNaseData[,motif_LF]),
    motif_RF = rowSums(DNaseData[,motif_RF]),

    flank_RF1 = rowSums(DNaseData[,flank_RF1]),
    flank_RF2 = rowSums(DNaseData[,flank_RF2]),
    flank_RF3 = rowSums(DNaseData[,flank_RF3]),
    flank_RF4 = rowSums(DNaseData[,flank_RF4]),
    flank_RF5 = rowSums(DNaseData[,flank_RF5]),

    flank_LR5 = rowSums(DNaseData[,flank_LR5]),
    flank_LR4 = rowSums(DNaseData[,flank_LR4]),
    flank_LR3 = rowSums(DNaseData[,flank_LR3]),
    flank_LR2 = rowSums(DNaseData[,flank_LR2]),
    flank_LR1 = rowSums(DNaseData[,flank_LR1]),

    motif_LR = rowSums(DNaseData[,motif_LR]),
    motif_RR = rowSums(DNaseData[,motif_RR]),

    flank_RR1 = rowSums(DNaseData[,flank_RR1]),
    flank_RR2 = rowSums(DNaseData[,flank_RR2]),
    flank_RR3 = rowSums(DNaseData[,flank_RR3]),
    flank_RR4 = rowSums(DNaseData[,flank_RR4]),
    flank_RR5 = rowSums(DNaseData[,flank_RR5]))

  M10 = data.frame(
    left2_fwd_sum = rowSums(DNaseData[,c(flank_LF5, flank_LF4)]),
    left1_fwd_sum = rowSums(DNaseData[,c(flank_LF3, flank_LF2, flank_LF1)]),
    motif_fwd_sum = rowSums(DNaseData[,c(motif_LF, motif_RF)]),
    right1_fwd_sum = rowSums(DNaseData[,c(flank_RF1, flank_RF2, flank_LF3)]),
    right2_fwd_sum = rowSums(DNaseData[,c(flank_RF4, flank_RF5)]),

    left2_rev_sum = rowSums(DNaseData[,c(flank_LR5, flank_LR4)]),
    left1_rev_sum = rowSums(DNaseData[,c(flank_LR3, flank_LR2, flank_LR1)]),
    motif_rev_sum = rowSums(DNaseData[,c(motif_LR, motif_RR)]),
    right1_rev_sum = rowSums(DNaseData[,c(flank_RR1, flank_RR2, flank_LR3)]),
    right2_rev_sum = rowSums(DNaseData[,c(flank_RR4, flank_RR5)]))

  if (combine_strands == 'motif') {
    M12 = M24[,c(1:12)] + M24[,c(24:13)]
  } else {
    # vertical combine
    M12 = M24[,c(1:12)] + M24[,c(13:24)]
  }
  colnames(M12) <- c('L5', 'L4', 'L3', 'L2', 'L1', 'motif_L', 'motif_R', 'R1', 'R2', 'R3', 'R4', 'R5')

  left2_sum = rowSums(M12[,c(1,2)])
  left1_sum = rowSums(M12[,c(3,4,5)])
  motif_sum = rowSums(M12[,c(6, 7)])
  right1_sum = rowSums(M12[,c(8,9,10)])
  right2_sum = rowSums(M12[,c(11,12)])

  M11 = data.frame(M12[,c(1:5)], motif_sum, M12[,c(8:12)])

  M5 = data.frame(left2_sum, left1_sum, motif_sum, right1_sum, right2_sum)

  M3 = data.frame(left1_sum, motif_sum, right1_sum)

  M2 = data.frame(left1_sum, right1_sum)

  M1 = left1_sum + right1_sum

  cuts_sum.l = list(M24 = M24,
                    M12 = M12,
                    M11 = M11,
                    M10 = M10,
                    M5 = M5,
                    M3 = M3,
                    M2 = M2,
                    M1 = M1)


  return(cuts_sum.l)

}

#' @title Normalize, binning and tranform DNase counts
#'
#' @description Normalize DNase tagcounts by library size,
#' binning DNase using MILLIPEDE binning and then take log2 or asinh transform
#' @param dnase_tagcounts DNase tagcount matrix
#' @param file_idxstats idxstats file generated by samtools
#' @param type_bin Type of MILLIPEDE binning scheme (Default: 'M5').
#' @param library_reference DNase-seq reference library size (Default: 100 million)
#' @param transform log2 or asinh transform
#'
#' @return
normalize_bin_dnase <- function(dnase_tagcounts, file_idxstats, type_bin='M5',
                                library_reference=1e8, transform='log2') {

  ## Count total mapped DNase reads (chr1:22)
  total_readsMapped <- readstats_total(file_idxstats, select_chr = TRUE)

  ## Normalize (scale) DNase read counts
  cat('Normalize DNase reads to', library_reference / 1e6, 'million \n')
  scaling_factor <- library_reference / total_readsMapped
  dnase_tagcounts_normalized = dnase_tagcounts * scaling_factor

  ## Bin DNase counts
  cat('Bin DNase data with', type_bin, 'binning \n')
  dnase_bins_normalized.df <- as.data.frame(sum_cuts_strand(dnase_tagcounts_normalized)[[type_bin]])

  if (transform == 'asinh') {
    cat(transform, 'transform \n')
    dnase_bins_normalized.df <- asinh(dnase_bins_normalized.df)
  } else if (transform == 'log2') {
    cat(transform, 'transform \n')
    dnase_bins_normalized.df <- log2(dnase_bins_normalized.df+1)
  }

  return(dnase_bins_normalized.df)

}

#' @title Normalize, binning and tranform DNase counts
#'
#' @param chip_counts ChIP-seq counts
#' @param file_idxstats idxstats file generated by samtools
#' @param library_reference ChIP-Seq reference library size (Default: 10 million)
#' @param transform log2 or asinh transform
#'
#' @return
normalize_chip <- function(chip_counts, file_idxstats, library_reference=1e7, transform='asinh') {

  ## Count total mapped DNase reads (chr1:22)
  total_readsMapped <- readstats_total(file_idxstats, select_chr = TRUE)

  ## Normalize (scale) ChIP-seq read counts
  cat('Normalize ChIP-seq reads to', library_reference / 1e6, 'million \n')
  scaling_factor <- library_reference / total_readsMapped
  chip_normalized <- chip_counts * scaling_factor

  if (transform == 'asinh') {
    cat(transform, 'transform \n')
    chip_normalized <- asinh(chip_normalized)
  } else if (transform == 'log2') {
    cat(transform, 'transform \n')
    chip_normalized <- log2(chip_normalized+1)
  }

  return(chip_normalized)

}

#' @title Get candidate sites by expanding FIMO motifs with flanking regions
#'
#' @param file_fimo Filename of FIMO results.
#' @param flank Size (bp) of flanking region on each side of motif. Default: 100.
#'
#' @return
#' @export
flank_fimo_sites <- function(file_fimo, flank=100) {

  cat('Load FIMO result and expand FIMO motifs with', flank, 'bp flanking regions. \n')

  fimo.df <- read.table(file_fimo, header = T, stringsAsFactors = F, comment.char = '', sep = '\t')

  if(nrow(fimo.df) < 1){
    cat('ERROR: No motif matches in FIMO result.\n', file=stderr())
    cleanup()
    q()
  }

  ## order sites
  chr_order <- paste0('chr', c(1:22, 'X','Y','M'))
  fimo.df <- fimo.df[fimo.df$sequence_name %in% chr_order,]
  fimo.df$sequence_name <- factor(fimo.df$sequence_name, chr_order, ordered=TRUE)
  fimo.df <- fimo.df[order(fimo.df$sequence_name, fimo.df$start, fimo.df$stop),]

  ## Collect fimo output's coordinates, pwm score, strand, p-value, q-value
  sites.df <- data.frame(chr = fimo.df$sequence_name,
                         start = fimo.df$start,
                         stop = fimo.df$stop,
                         name = paste0('site', c(1:nrow(fimo.df))),
                         # name = paste0(fimo.df$sequence_name, ':', fimo.df$start, '-', fimo.df$end, '.', fimo.df$strand),
                         score = fimo.df$score,
                         strand = fimo.df$strand,
                         p.value = fimo.df$p.value)

  ## FIMO output are 1-based, convert to BED format [Start: end) 0-based coordinates
  sites.df$start <- sites.df$start -1
  sites.df$stop <- sites.df$stop

  ## Expand flanking regions
  sites.df$start <- sites.df$start - flank
  sites.df$stop <- sites.df$stop + flank

  ## Filter out sites with start positions < 0 after flanking
  sites.df <- sites.df[sites.df$start >= 0, ]

  return(sites.df)

}

#' @title Get candidate sites by expanding FIMO motifs with flanking regions
#'
#' @param sites data frame of the candidate sites (BED format).
#' @param file_mapability Filename of the mapability reference file in bigWig format.
#' @param path_bigWigAverageOverBed Path to bigWigAverageOverBed executable.
#'
#' @return
#' @export
compute_mapability <- function(sites,
                               file_mapability,
                               path_bigWigAverageOverBed='bigWigAverageOverBed') {

  if ( system( paste(path_bigWigAverageOverBed,'--help') , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( 'ERROR: bigWigAverageOverBed could not be executed, set path_bigWigAverageOverBed\n' , sep='', file=stderr() )
    cleanup()
    q()
  }

  if( !file.exists(file_mapability) ){
    cat( 'Cannot find the mapability bigWig file, set file_mapability\n' , sep='', file=stderr() )
    cleanup()
    q()
  }

  sites_tmp <- sites[,1:6]
  sites_tmp$name <- paste('site', c(1:nrow(sites_tmp)), sep = '')
  sites_tmp$score <- 0

  tmp_file_sites <- tempfile('sites')
  write.table(sites_tmp, tmp_file_sites, sep = '\t', quote = F, row.names = F, col.names = F)

  cat('Compute mapability ... \n')
  tmp_file_mapability <- tempfile('mapability')
  cmd <- paste(path_bigWigAverageOverBed, file_mapability, tmp_file_sites, tmp_file_mapability)
  system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

  mapability <- read.table(tmp_file_mapability)[,5]

  file.remove(c(tmp_file_sites, tmp_file_mapability))

  return(mapability)

}

#' @title Filter out sites in blacklist regions
#'
#' @param sites data frame of the candidate sites (BED format).
#' @param file_blacklist Filename of the blacklist regions
#' @param path_bedtools Path to bedtools executable.
#'
#' @return
#' @export
filter_blacklist <- function(sites,
                             file_blacklist,
                             path_bedtools='bedtools') {

  if ( system( paste(path_bedtools,'--help') , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
    cat( 'ERROR: bedtools could not be executed, set path_bedtools\n' , sep='', file=stderr() )
    cleanup()
    q()
  }

  if( !file.exists(file_blacklist) ){
    cat( 'Cannot find the blacklist file, check file_blacklist\n' , sep='', file=stderr() )
    cleanup()
    q()
  }

  sites_tmp <- sites[,1:6]
  sites_tmp$name <- paste('site', c(1:nrow(sites_tmp)), sep = '')
  sites_tmp$score <- 0

  tmp_file_sites <- tempfile('sites')
  write.table(sites_tmp, tmp_file_sites, sep = '\t', quote = F, row.names = F, col.names = F)

  cat('Filter blacklist regions \n')
  tmp_file_sites_filtered <- tempfile('sites_filtered')
  cmd <- paste(path_bedtools, 'intersect -v -a', tmp_file_sites, '-b', file_blacklist, '>', tmp_file_sites_filtered)
  system( cmd, ignore.stdout=SYS_PRINT, ignore.stderr=SYS_PRINT )

  sites_filtered <- read.table(tmp_file_sites_filtered)
  colnames(sites_filtered) <- colnames(sites)

  file.remove(c(tmp_file_sites, tmp_file_sites_filtered))

  return(sites_filtered)

}


#' @title Process candidate sites from FIMO result
#'
#' @param file_fimo Filename of FIMO result.
#' @param file_sites Filenmae of candidate sites.
#' @param flank Size (bp) of flanking region on each side of motif. Default: 100.
#' @param thresh_pValue FIMO p-value threshold.
#' @param compute_mapability If TRUE, compute mapability for candidate sites.
#' @param file_mapability Filename of the mapability reference file in bigWig format.
#' @param filter_blacklist If TRUE, filter out sites in blacklist regions.
#' @param file_blacklist Filename of the blacklist regions
#' @param path_bigWigAverageOverBed Path to bigWigAverageOverBed executable.
#' @param path_bedtools Path to bedtools executable.
#'
#' @export
process_candidate_sites <- function(file_fimo,
                                    file_sites,
                                    flank=100,
                                    thresh_pValue=1e-5,
                                    compute_mapability=FALSE,
                                    file_mapability=NA,
                                    filter_blacklist=FALSE,
                                    file_blacklist=NA,
                                    path_bigWigAverageOverBed='bigWigAverageOverBed',
                                    path_bedtools='bedtools') {

  options(scipen=999) ## suppress scientific notations

  ## Get candidate sites by expanding FIMO motifs with flanking regions
  sites.df <- flank_fimo_sites(file_fimo, flank)

  ## Select candidate sites with FIMO p-value < thresh_pValue
  sites.df <- sites.df[which(as.numeric(sites.df$p.value) < as.numeric(thresh_pValue)), ]
  cat(nrow(sites.df), 'sites with p-value <', thresh_pValue, '\n')

  if(compute_mapability) {
    ## Compute mapability for sites
    sites.df$mapability <- compute_mapability(sites.df, file_mapability, path_bigWigAverageOverBed)
  }

  if(filter_blacklist) {
    sites.df <- filter_blacklist(sites.df, file_blacklist, path_bedtools)
  }

  if(!dir.exists(dirname(file_sites))){
    dir.create(dirname(file_sites), showWarnings = F, recursive = T)
  }

  colnames(sites.df)[1] <- paste0('#', colnames(sites.df)[1])
  write.table(sites.df, file_sites, sep = '\t', quote = F, row.names = F, col.names = T)
  cat(nrow(sites.df), 'candidate sites written in', filename_sites, '\n')

}
