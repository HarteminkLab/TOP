
#' Get candidate sites using FIMO motifs with flanking regions
#'
#' @param fimo_file FIMO result .txt file
#' @param flank Flanking region (bp) around motif matches (default: 100)
flank_fimo_sites <- function(fimo_file, flank=100) {

  if( !file.exists(fimo_file) || file.size(fimo_file) == 0 ){
    stop(paste(fimo_file, 'file does not exist or is empty. Probably no motif matches were found.'))
  }

  cat('Load FIMO result and flank motifs by', flank, 'bp... \n')
  fimo.df <- read.table(fimo_file, header = TRUE, stringsAsFactors = FALSE,
                        comment.char = '', sep = '\t')
  # Sort sites
  chr_order <- paste0('chr', c(1:22, 'X','Y','M'))
  fimo.df <- fimo.df[fimo.df$sequence_name %in% chr_order,]
  fimo.df$sequence_name <- factor(fimo.df$sequence_name, chr_order, ordered=TRUE)
  fimo.df <- fimo.df[order(fimo.df$sequence_name, fimo.df$start, fimo.df$stop),]

  ## Prepare candidate sites
  sites.df <- data.frame(chr = fimo.df$sequence_name,
                         start = fimo.df$start,
                         stop = fimo.df$stop,
                         name = paste0('site', c(1:nrow(fimo.df))),
                         score = fimo.df$score,
                         strand = fimo.df$strand,
                         p.value = fimo.df$p.value)

  # FIMO output are 1-based, convert to BED format 0-based coordinates
  sites.df$start <- sites.df$start -1
  sites.df$stop <- sites.df$stop

  # Expand flanking regions
  sites.df$start <- sites.df$start - flank
  sites.df$stop <- sites.df$stop + flank

  # Filter out sites with start positions < 0 after flanking
  sites.df <- sites.df[sites.df$start >= 0, ]

  return(sites.df)

}

#' Filter sites in ENCODE blacklist regions
#'
#' @param sites.df A data frame of candidate binding sites
#' @param blacklist_file ENCODE blacklist file
#'
#' @importFrom data.table fread
#' @import GenomicRanges
filter_blacklist <- function(sites.df, blacklist_file) {

  if( !file.exists(blacklist_file) ){
    stop( 'Cannot location the blacklist file!' )
  }

  blacklist <- fread(blacklist_file)
  colnames(blacklist)[1:3] <- c('chr', 'start', 'end')
  blacklist.gr <- makeGRangesFromDataFrame(blacklist)

  # Filter sites overlapping with blacklist regions
  sites.gr <- makeGRangesFromDataFrame(sites.df, keep.extra.columns = T)
  in.blacklist <- which(countOverlaps(query = sites.gr, subject = blacklist.gr))
  cat('Filter out', length(in.blacklist), 'sites overlapping with blacklist regions. \n')
  sites.df <- sites.df[-in.blacklist, ]

  return(sites.df)

}

# Compute mapability for candidate sites
#' Title
#'
#' @param sites A data frame of candidate binding sites with the first
#' 6 columns the same as in the BED format
#' @param mapability_file ENCODE mapability bigWig file
#' @param bigWigAverageOverBed_path path of bigWigAverageOverBed executable
#' @importFrom data.table fread fwrite
#'
compute_mapability <- function(sites,
                               mapability_file,
                               bigWigAverageOverBed_path='bigWigAverageOverBed') {

  if( !file.exists(mapability_file) ){
    stop( 'Cannot locate the mapability bigWig file!' )
  }

  if ( system(paste(bigWigAverageOverBed_path,'--help'), ignore.stdout=T, ignore.stderr=T ) != 0 ) {
    stop( 'bigWigAverageOverBed could not be executed, set bigWigAverageOverBed_path!' )
  }

  sites_tmp <- sites[,1:6]
  sites_tmp$name <- paste('site', c(1:nrow(sites_tmp)), sep = '')
  sites_tmp$score <- 0

  tmp_mapability_filesites <- tempfile('sites')
  fwrite(sites_tmp, tmp_mapability_filesites, sep = '\t')

  cat('Compute mapability ... \n')
  tmp_mapability_file <- tempfile('mapability')
  cmd <- paste(bigWigAverageOverBed_path,
               mapability_file, tmp_mapability_filesites, tmp_mapability_file)
  system(cmd)

  mapability <- fread(tmp_mapability_file)[,5]

  file.remove(c(tmp_mapability_filesites, tmp_mapability_file))

  return(mapability)

}



#' @title Full process to get candidate sites from FIMO result
#'
#' @param fimo_file Filename of FIMO result.
#' @param flank Flanking region (bp) around motif matches (default: 100)
#' @param thresh_pValue FIMO p-value threshold.
#' @param blacklist_file Filename of the blacklist regions
#' @param mapability_file Filename of the mapability reference file in bigWig format.
#' @param out_file Filename of processed candidate sites.
#' @param bigWigAverageOverBed_path Path to bigWigAverageOverBed executable.
#' This is only needed for computing mapability.
#' @importFrom data.table fread fwrite
#'
process_candidate_sites <- function(fimo_file,
                                    flank=100,
                                    thresh_pValue=1e-5,
                                    blacklist_file,
                                    mapability_file,
                                    out_file,
                                    bigWigAverageOverBed_path='bigWigAverageOverBed') {

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

  if(is.null(out_file)){
    if(!dir.exists(dirname(out_file))){
      dir.create(dirname(out_file), showWarnings = F, recursive = T)
    }
    colnames(sites.df)[1] <- paste0('#', colnames(sites.df)[1])
    fwrite(sites.df, out_file, sep = '\t')
    cat('Candidate sites output to', out_file, '\n')
  }

  return(sites.df)
}
