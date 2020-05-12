#' Match DNase tagcount matrices for candidate sites
#'
#' @param file_sites filename for candidate sites
#' @param file_dnase_tagcount_fwd filename for DNase tagcount in forward strand
#' @param file_dnase_tagcount_rev filename for DNase tagcount in reverse strand
#' @param file_dnase_matrix_fwd filename for DNase tagcount matrix in forward strand
#' @param file_dnase_matrix_rev filename for DNase tagcount matrix in reverse strand
#'
#' @export
match_tagcount_sites <- function(file_sites,
                                 file_dnase_tagcount_fwd, file_dnase_tagcount_rev,
                                 file_dnase_matrix_fwd, file_dnase_matrix_rev) {

  cat('Match matrix of', file_dnase_matrix_fwd, '... \n')
  cmd_fwd <- paste('cut -f 1-4', file_sites, '| bwtool extract bed stdin',
                   file_dnase_tagcount_fwd, file_dnase_matrix_fwd, '-fill=0 -decimals=0 -tabs')
  system(cmd_fwd)

  cat('Match matrix of', file_dnase_matrix_rev, '... \n')
  cmd_rev <- paste('cut -f 1-4', file_sites, '| bwtool extract bed stdin',
                   file_dnase_tagcount_rev, file_dnase_matrix_rev, '-fill=0 -decimals=0 -tabs')
  system(cmd_rev)

  # Flip the tagcounts generated from bwtool for motifs on the reverse strand
  rev_tagcount_bwtool(file_sites, file_dnase_matrix_fwd, file_dnase_matrix_rev)
}

#' Flip the tagcounts generated from bwtool for motifs on the reverse (minus) strand
#'
#' @param file_sites filename for candidate sites
#' @param file_dnase_matrix_fwd filename for DNase tagcount matrix in forward strand
#' @param file_dnase_matrix_rev filename for DNase tagcount matrix in reverse strand
#'
#' @return
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
