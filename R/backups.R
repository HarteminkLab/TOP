
# Counts DNase-seq or ATAC-seq cuts along the genome (without using bedtools)
# This is an alternative function to the count_genome_cuts(),
# but without using \code{bedtools}.
# This is often slower than count_genome_cuts() and may require more memory
# for large BAM files.
count_genome_cuts_nobedtools <- function(bam_file,
                                         chrom_size_file,
                                         data_type=c('DNase', 'ATAC'),
                                         shift_ATAC = FALSE,
                                         shift_ATAC_bases=c(4L,-4L),
                                         outdir = dirname(bam_file),
                                         outname,
                                         bedGraphToBigWig_path='bedGraphToBigWig'){

  # Checking input arguments
  if ( Sys.which(bedGraphToBigWig_path) == '' )
    stop( 'bedGraphToBigWig could not be executed. Please install bedGraphToBigWig and set bedGraphToBigWig_path.' )

  data_type <- match.arg(data_type)

  coverage.gr <- read_bam_cuts(bam_file, data_type, shift_ATAC, shift_ATAC_bases, return_type = 'coverage')
  pos_coverage.gr <- coverage.gr$pos
  neg_coverage.gr <- coverage.gr$neg

  rm(coverage.gr)

  if(!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  if(missing(outname))
    outname <- tools::file_path_sans_ext(basename(bam_file))

  # cuts coverage on + strand counts
  cat('Processing genome cuts on + strand ...\n')
  pos_coverage.df <- data.frame(chr = as.character(seqnames(pos_coverage.gr)),
                                start = start(pos_coverage.gr) - 1,
                                end = end(pos_coverage.gr),
                                count = pos_coverage.gr$score)
  pos_coverage.df <- pos_coverage.df[pos_coverage.df$count > 0, ]
  pos_coverage.df <- pos_coverage.df[with(pos_coverage.df, order(chr, start)), ]
  bedgraph_file <- file.path(outdir, paste0(outname, '.fwd.genomecounts.bedGraph'))
  data.table::fwrite(pos_coverage.df, bedgraph_file,
                     sep='\t', col.names = FALSE, scipen = 999)
  bedGraphToBigWig(bedgraph_file, chrom_size_file, bedGraphToBigWig_path)
  unlink(bedgraph_file)

  # cuts coverage on - strand counts
  cat('Processing genome cuts on - strand ...\n')
  neg_coverage.df <- data.frame(chr = as.character(seqnames(neg_coverage.gr)),
                                start = start(neg_coverage.gr) - 1,
                                end = end(neg_coverage.gr),
                                count = neg_coverage.gr$score)
  neg_coverage.df <- neg_coverage.df[neg_coverage.df$count > 0, ]
  neg_coverage.df <- neg_coverage.df[with(neg_coverage.df, order(chr, start)), ]
  bedgraph_file <- file.path(outdir, paste0(outname, '.rev.genomecounts.bedGraph'))
  data.table::fwrite(neg_coverage.df, bedgraph_file,
                     sep='\t', col.names = FALSE, scipen = 999)
  bedGraphToBigWig(bedgraph_file, chrom_size_file, bedGraphToBigWig_path)
  unlink(bedgraph_file)

}


# Reads ATAC cuts for ATAC-seq paired-end reads
read_bam_cuts_ATACreadpairs <- function(bam_file,
                                        select_NFR_fragments = FALSE,
                                        shift_ATAC = FALSE,
                                        shift_ATAC_bases=c(4L,-4L),
                                        return_type = c('cuts', 'coverage')){

  return_type <- match.arg(return_type)
  cat('Reading genome cuts for', bam_file, '...\n')

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

  # Extract 5' end position and shift reads
  if(shift_ATAC){
    cat(sprintf('Shifting ATAC-seq reads by %d and %d bp...\n',
                shift_ATAC_bases[1], shift_ATAC_bases[2]))

    read1_5ends <- resize(granges(read1), fix = 'start', 1)
    read1_pos_cuts <- shift(granges(read1_5ends[strand(read1_5ends) == '+']), shift_ATAC_bases[1])
    read1_neg_cuts <- shift(granges(read1_5ends[strand(read1_5ends) == '-']), shift_ATAC_bases[2])

    read2_5ends <- resize(granges(read2), fix = 'start', 1)
    read2_pos_cuts <- shift(granges(read2_5ends[strand(read2_5ends) == '+']), shift_ATAC_bases[1])
    read2_neg_cuts <- shift(granges(read2_5ends[strand(read2_5ends) == '-']), shift_ATAC_bases[2])

    cuts <- c(read1_pos_cuts, read1_neg_cuts,
              read2_pos_cuts, read2_neg_cuts)
  }else{
    read1_5ends <- resize(granges(read1), fix = 'start', 1)
    read2_5ends <- resize(granges(read2), fix = 'start', 1)
    cuts <- c(read1_5ends, read2_5ends)
  }

  if(return_type == 'coverage'){
    cat('Counting genome cuts coverage .. \n')
    pos_coverage.gr <- GRanges(coverage(cuts[strand(cuts) == '+']), strand = '+')
    neg_coverage.gr <- GRanges(coverage(cuts[strand(cuts) == '-']), strand = '-')
    cuts_coverage.gr <- GRangesList(pos = pos_coverage.gr,
                                    neg = neg_coverage.gr)
    return(cuts_coverage.gr)
  }else if(return_type == 'cuts'){
    return(cuts)
  }

}


# Reads DNase or ATAC cuts or coverage of cuts from BAM file
read_bam_cuts <- function(bam_file,
                          data_type=c('DNase', 'ATAC'),
                          shift_ATAC = FALSE,
                          shift_ATAC_bases=c(4L,-4L),
                          return_type = c('cuts', 'coverage')){

  data_type <- match.arg(data_type)
  return_type <- match.arg(return_type)
  cat('Reading genome cuts for', bam_file, '...\n')

  # Read ATAC-seq alignments
  reads <- GenomicAlignments::readGAlignments(bam_file)

  # Extract 5' end position and shift based on strand
  if (data_type == 'ATAC' && shift_ATAC == TRUE) {
    cat('Shifting ATAC-seq reads ...\n')
    cuts <- GenomicRanges::resize(GenomicRanges::granges(reads), fix = 'start', 1)
    pos_cuts <- GenomicRanges::shift(GenomicRanges::granges(cuts[strand(cuts) == '+']), shift_ATAC_bases[1])
    neg_cuts <- GenomicRanges::shift(GenomicRanges::granges(cuts[strand(cuts) == '-']), shift_ATAC_bases[2])
    cuts <- c(pos_cuts, neg_cuts)
  }else{
    cuts <- GenomicRanges::resize(GenomicRanges::granges(reads), fix = 'start', 1)
  }

  if(return_type == 'coverage'){
    cat('Counting genome cuts coverage .. \n')
    pos_coverage.gr <- GRanges(coverage(cuts[strand(cuts) == '+']), strand = '+')
    neg_coverage.gr <- GRanges(coverage(cuts[strand(cuts) == '-']), strand = '-')
    cuts_coverage.gr <- GRangesList(pos = pos_coverage.gr,
                                    neg = neg_coverage.gr)
    return(cuts_coverage.gr)
  }else if(return_type == 'cuts'){
    return(cuts)
  }

}


