
#' Get candidate sites using FIMO motifs with flanking regions
#'
#' @param fimo_file FIMO result .txt file
#' @param flank Flanking region (bp) around motif matches (default: 100)
#' @importFrom data.table fread
#' @export
#'
flank_fimo_sites <- function(fimo_file, flank=100) {

  if( !file.exists(fimo_file) ){
    stop(paste(fimo_file, 'file does not exist or is empty!'))
  }

  cat('Load FIMO result... \n')
  fimo.df <- as.data.frame(fread(fimo_file, sep ='\t'))

  # Sort sites
  chr_order <- paste0('chr', c(1:22, 'X','Y','M'))
  fimo.df <- fimo.df[fimo.df$sequence_name %in% chr_order,]
  fimo.df$sequence_name <- factor(fimo.df$sequence_name, chr_order, ordered=TRUE)
  fimo.df <- fimo.df[order(fimo.df$sequence_name, fimo.df$start, fimo.df$stop),]

  cat('Flank motif matches by', flank, 'bp... \n')
  ## Prepare candidate sites
  sites.df <- data.frame(chr = fimo.df$sequence_name,
                         start = fimo.df$start,
                         end = fimo.df$stop,
                         name = paste0('site', c(1:nrow(fimo.df))),
                         pwm.score = fimo.df$score,
                         strand = fimo.df$strand,
                         p.value = fimo.df$`p-value`)

  # FIMO output are 1-based, convert to BED format 0-based coordinates
  sites.df$start <- sites.df$start -1
  sites.df$end <- sites.df$end

  # Expand flanking regions
  sites.df$start <- sites.df$start - flank
  sites.df$end <- sites.df$end + flank

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
#' @export
#'
filter_blacklist <- function(sites.df, blacklist_file) {

  if( !file.exists(blacklist_file) ){
    stop( 'Cannot location the blacklist file!' )
  }

  blacklist <- fread(blacklist_file)
  colnames(blacklist)[1:3] <- c('chr', 'start', 'end')
  blacklist.gr <- makeGRangesFromDataFrame(blacklist)

  # Filter sites overlapping with blacklist regions
  sites.gr <- makeGRangesFromDataFrame(sites.df, keep.extra.columns = T)
  in.blacklist <- which(countOverlaps(query = sites.gr, subject = blacklist.gr)>0)
  cat('Filter out', length(in.blacklist), 'sites in blacklist regions. \n')
  sites.df <- sites.df[-in.blacklist, ]

  return(sites.df)

}

# Filter mapability for candidate sites
#' Title
#'
#' @param sites.df A data frame of candidate binding sites with the first
#' 6 columns the same as in the BED format
#' @param mapability_file ENCODE mapability bigWig file
#' @param thresh_mapability Mpability filter threshold
#' @param bigWigAverageOverBed_path path of bigWigAverageOverBed executable
#' @importFrom data.table fread fwrite
#'
#' @export
#'
filter_mapability <- function(sites.df=NULL,
                              mapability_file=NULL,
                              thresh_mapability=0.8,
                              bigWigAverageOverBed_path='bigWigAverageOverBed') {

  if( !file.exists(mapability_file) ){
    stop( 'Cannot locate the mapability bigWig file!' )
  }

  if ( system(paste(bigWigAverageOverBed_path,'--help'), ignore.stdout=T, ignore.stderr=T ) != 0 ) {
    stop( 'bigWigAverageOverBed could not be executed, set bigWigAverageOverBed_path!' )
  }

  sites_tmp.df <- sites.df[,1:6]
  sites_tmp.df$name <- paste('site', c(1:nrow(sites_tmp.df)), sep = '')
  sites_tmp.df$pwm.score <- 0

  tmp_mapability_filesites <- tempfile('sites')
  fwrite(sites_tmp.df, tmp_mapability_filesites, sep = '\t')

  cat('Compute mapability ... \n')
  tmp_mapability_file <- tempfile('mapability')
  system(paste(bigWigAverageOverBed_path, mapability_file, tmp_mapability_filesites, tmp_mapability_file))

  sites.df$mapability <- fread(tmp_mapability_file)[,5]

  sites.df <- sites.df[sites.df$mapability > thresh_mapability, ]

  file.remove(c(tmp_mapability_filesites, tmp_mapability_file))

  return(sites.df)

}



#' @title Full process to get candidate sites from FIMO result
#'
#' @param fimo_file Filename of FIMO result.
#' @param flank Flanking region (bp) around motif matches (default: 100)
#' @param thresh_pValue FIMO p-value threshold.
#' @param thresh_pwmscore FIMO PWM score threshold.
#' @param blacklist_file Filename of the blacklist regions
#' @param mapability_file Filename of the mapability reference file in bigWig format.
#' @param thresh_mapability Mapability threshold (default: 0.8,
#' include sites map-able at least 80% positions).
#' @param out_file Filename of processed candidate sites.
#' @param bigWigAverageOverBed_path Path to bigWigAverageOverBed executable.
#' This is only needed for computing mapability.
#' @importFrom data.table fread fwrite
#'
#' @export
#'
process_candidate_sites <- function(fimo_file=NULL,
                                    flank=100,
                                    thresh_pValue=1e-5,
                                    thresh_pwmscore=0,
                                    blacklist_file=NULL,
                                    mapability_file=NULL,
                                    thresh_mapability=0.8,
                                    out_file=NULL,
                                    bigWigAverageOverBed_path='bigWigAverageOverBed') {

  # Get candidate sites from FIMO motif matches and add flanking regions
  sites.df <- flank_fimo_sites(fimo_file, flank)

  # Filter candidate sites by FIMO p-value
  sites.df <- sites.df[which(as.numeric(sites.df$p.value) < as.numeric(thresh_pValue)), ]
  cat('Select candidate sites with FIMO p-value <', thresh_pValue, '\n')

  # Filter candidate sites by FIMO PWM score
  sites.df <- sites.df[which(as.numeric(sites.df$pwm.score) > as.numeric(thresh_pwmscore)), ]
  cat('Select candidate sites with PWM score >', thresh_pwmscore, '\n')

  # Filter candidate sites in ENCODE blacklist
  if(!is.null(blacklist_file)) {
    sites.df <- filter_blacklist(sites.df, blacklist_file)
  }

  # Filter candidate sites by mapability
  if(!is.null(mapability_file)) {
    sites.df <- filter_mapability(sites.df, mapability_file,
                                  thresh_mapability, bigWigAverageOverBed_path)
  }

  if(!is.null(out_file)){
    if(!dir.exists(dirname(out_file))){
      dir.create(dirname(out_file), showWarnings = F, recursive = T)
    }
    tmp.df <- sites.df
    colnames(tmp.df)[1] <- paste0('#', colnames(tmp.df)[1])
    fwrite(tmp.df, out_file, sep = '\t')
    cat('Save candidate sites at', out_file, '\n')
  }

  return(sites.df)
}
