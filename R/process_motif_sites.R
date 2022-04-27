
#' @title Runs \code{FIMO} to scan for motif matches
#' @description Runs \code{FIMO} to scan for motif matches along the genome.
#' @param motif_file Motif file in MEME format.
#' @param sequence_file Genome sequence file in FASTA format.
#' @param outname Output file name.
#' @param outdir Output directory.
#' @param thresh_pValue \code{FIMO} option \code{thresh} for p-value threshold.
#' @param background Option for background model:
#' \sQuote{default}: uses \code{FIMO} default background setting;
#' \sQuote{motif}: uses the 0-order letter frequencies contained in the motif file;
#' \sQuote{uniform}: uses uniform letter frequencies;
#' \sQuote{file}: uses the file specified in \code{background_file}.
#' @param background_file Path to a file in Markov Background Model Format.
#' @param skip_matched_sequence \code{FIMO} option \code{skip_matched_sequence}.
#' Turns off output of the sequence of motif matches.
#' This speeds up processing considerably.
#' @param max_strand \code{FIMO} option \code{max_strand}.
#' If matches on both strands at a
#' given position satisfy the output threshold,
#' only report the match for the strand with the higher score.
#' If the scores are tied, the matching strand is chosen at random.
#' @param max_stored_scores The maximum number of stored matches.
#' @param options Other options for \code{FIMO}.
#' @param verbosity A number of the verbosity level (from 1 to 5).
#' If set to 1 (quiet) then it will only output error messages, in contrast, the
#' other extreme 5 (dump) outputs lots of mostly useless information.
#' @param fimo_path Path to \code{fimo} command line executable.
#' @export
#' @examples
#' \dontrun{
#' fimo_motif_matches(motif_file='motifID.meme',
#'                    sequence_file='hg38.fa',
#'                    thresh_pValue=1e-5,
#'                    outname='motifID_1e-5.fimo.txt',
#'                    fimo_path='fimo')
#' }
fimo_motif_matches <- function(motif_file,
                               sequence_file,
                               outname='fimo.txt',
                               outdir=dirname(outname),
                               thresh_pValue=1e-5,
                               background=c('default', 'motif', 'uniform', 'file'),
                               background_file,
                               skip_matched_sequence=TRUE,
                               max_strand=FALSE,
                               max_stored_scores=100000,
                               options='',
                               verbosity=2,
                               fimo_path='fimo') {

  background <- match.arg(background)

  if ( Sys.which(fimo_path) == '' ) {
    stop( 'fimo could not be executed. Please install fimo and set fimo_path.' )
  }

  if(max_strand){
    max_strand <- '--max-strand'
  }else{
    max_strand <- ''
  }

  if(skip_matched_sequence){
    skip_matched_sequence <- '--skip-matched-sequence'
  }else{
    skip_matched_sequence <- ''
  }

  if(background == 'default'){
    bfile <- ''
  }else if(background == 'motif'){
    bfile <- '--bfile --motif--'
  }else if(background == 'uniform'){
    bfile <- '--bfile --uniform--'
  }else if(background == 'file'){
    if(missing(background_file)){
      stop('Please specify the path to the background file!')
    }else if(!file.exists(background_file)){
      stop(paste('Background file', background_file, 'cannot be found!'))
    }else{
      bfile <- paste('--bfile', background_file)
    }
  }

  if(missing(max_stored_scores)){
    max_stored_scores <- 100000
  }

  if(!verbosity %in% c(1,2,3,4,5)){
    stop('verbosity must be 1, 2, 3, 4, or 5!')
  }

  if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }

  cmd <- paste('fimo --text',
               bfile,
               '--thresh', thresh_pValue,
               skip_matched_sequence,
               max_strand,
               '--oc', outdir,
               '--verbosity', verbosity,
               '--max-stored-scores', max_stored_scores,
               options,
               motif_file, sequence_file,
               '>', outname)
  cat(paste0('Run FIMO command: \n', cmd, '\n'))

  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

}


#' @title Obtains and filters candidate sites from FIMO result
#' @description Gets candidate sites from FIMO motif matches and
#' add flanking regions, and filters the candidate sites by different thresholds,
#' and filters out sites in blacklist regions.
#' @param fimo_file \code{FIMO} result \sQuote{.txt} file
#' @param flank Flanking region (bp) around motif matches (default: 100)
#' @param thresh_pValue \code{FIMO} p-value threshold (default: 1e-5)
#' @param thresh_pwmscore \code{FIMO} PWM score threshold (default: 0)
#' @param blacklist_file Filename of the blacklist regions (default: NULL)
#' @param mapability_file Filename of the mapability reference file
# in bigWig format (default: NULL)
#' @param thresh_mapability Mapability threshold (default: 0.8,
#' candidate sites need to be mapable at least 80% positions).
#' This is only used when the \code{mapability_file} is available.
#' @param bigWigAverageOverBed_path Path to bigWigAverageOverBed executable
#' (only needed when filtering mapability).
#' @return A data frame of processed candidate binding sites
#' passing the threshold filtering.
#' @export
#' @examples
#' \dontrun{
#' sites <- process_candidate_sites(fimo_file='fimo.txt',
#'                                  thresh_pValue=1e-5,
#'                                  blacklist_file='blacklist.hg38.bed.gz')
#' }
process_candidate_sites <- function(fimo_file,
                                    flank=100,
                                    thresh_pValue=1e-5,
                                    thresh_pwmscore=0,
                                    blacklist_file=NULL,
                                    mapability_file=NULL,
                                    thresh_mapability=0.8,
                                    bigWigAverageOverBed_path='bigWigAverageOverBed') {

  if( !file.exists(fimo_file) || file.size(fimo_file) == 0 ){
    stop(paste(fimo_file, 'file does not exist or is empty!'))
  }

  # Gets candidate sites from FIMO motif matches and add flanking regions
  sites <- flank_fimo_sites(fimo_file, flank)

  # Filters candidate sites by FIMO p-value
  sites <- sites[which(as.numeric(sites$p.value) < as.numeric(thresh_pValue)), ]
  cat('Select candidate sites with FIMO p-value <', thresh_pValue, '\n')

  # Filters candidate sites by PWM score
  sites <- sites[which(as.numeric(sites$pwm.score) > as.numeric(thresh_pwmscore)), ]
  cat('Select candidate sites with PWM score >', thresh_pwmscore, '\n')

  # Filters candidate sites in ENCODE blacklist
  if(!is.null(blacklist_file)) {
    sites <- filter_blacklist(sites, blacklist_file)
  }

  # Filters candidate sites by mapability
  if(!is.null(mapability_file)) {
    sites <- filter_mapability(sites, mapability_file,
                               thresh_mapability, bigWigAverageOverBed_path)
  }

  return(sites)
}

# Gets candidate sites using FIMO motifs with flanking regions
flank_fimo_sites <- function(fimo_file, flank=100) {

  if( !file.exists(fimo_file) || file.size(fimo_file) == 0 ){
    stop(paste(fimo_file, 'file does not exist or is empty!'))
  }

  cat('Load FIMO result... \n')
  fimo <- as.data.frame(fread(fimo_file, sep ='\t'))

  # Sorts sites
  chr_order <- paste0('chr', c(1:22, 'X','Y','M'))
  fimo <- fimo[fimo$sequence_name %in% chr_order,]
  fimo$sequence_name <- factor(fimo$sequence_name, chr_order, ordered=TRUE)
  fimo <- fimo[order(fimo$sequence_name, fimo$start, fimo$stop),]

  cat('Flank motif matches by', flank, 'bp... \n')
  sites <- data.frame(chr = fimo$sequence_name,
                      start = fimo$start,
                      end = fimo$stop,
                      name = paste0('site', c(1:nrow(fimo))),
                      pwm.score = fimo$score,
                      strand = fimo$strand,
                      p.value = fimo$`p-value`)

  # FIMO output are 1-based, convert to BED format 0-based coordinates
  sites$start <- sites$start -1
  sites$end <- sites$end

  # Expands flanking regions
  sites$start <- sites$start - flank
  sites$end <- sites$end + flank

  # Filters out sites with start positions < 0 after flanking
  sites <- sites[sites$start >= 0, ]

  return(sites)

}

# Filters candidate sites in blacklist regions
filter_blacklist <- function(sites, blacklist_file) {

  if( !file.exists(blacklist_file) ){
    stop( 'Cannot location the blacklist file!' )
  }

  blacklist <- fread(blacklist_file)
  colnames(blacklist)[1:3] <- c('chr', 'start', 'end')
  blacklist.gr <- makeGRangesFromDataFrame(blacklist)

  # Filter sites overlapping with blacklist regions
  sites.gr <- makeGRangesFromDataFrame(sites, keep.extra.columns = T)
  in.blacklist <- countOverlaps(query = sites.gr, subject = blacklist.gr)>0
  cat('Filter out', sum(in.blacklist), 'sites in blacklist regions. \n')
  if(any(in.blacklist)){
    sites <- sites[!in.blacklist, ]
  }

  return(sites)

}

# Filters candidate sites by mapability
filter_mapability <- function(sites,
                              mapability_file,
                              thresh_mapability=0.8,
                              bigWigAverageOverBed_path='bigWigAverageOverBed') {

  if( !file.exists(mapability_file) ){
    stop( 'Cannot locate the mapability bigWig file!' )
  }

  if ( system(paste(bigWigAverageOverBed_path,'--help'), ignore.stdout=T, ignore.stderr=T ) != 0 ) {
    stop( 'bigWigAverageOverBed could not be executed. Please install bigWigAverageOverBed and set bigWigAverageOverBed_path.' )
  }

  sites_tmp <- sites[,1:6]
  sites_tmp$name <- paste('site', c(1:nrow(sites_tmp)), sep = '')
  sites_tmp$pwm.score <- 0

  tmp_mapability_sites <- tempfile('sites')
  fwrite(sites_tmp, tmp_mapability_sites, sep = '\t', col.names = FALSE, scipen = 999)

  cat('Compute mapability ... \n')
  tmp_mapability_file <- tempfile('mapability')
  cmd <- paste(bigWigAverageOverBed_path, mapability_file, tmp_mapability_sites, tmp_mapability_file)
  if(.Platform$OS.type == 'windows') shell(cmd) else system(cmd)

  sites$mapability <- fread(tmp_mapability_file)[,5]

  sites <- sites[sites$mapability > thresh_mapability, ]

  unlink(c(tmp_mapability_sites, tmp_mapability_file))

  return(sites)

}


