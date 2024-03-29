% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_dnase_atac_data.R
\name{normalize_bin_transform_counts}
\alias{normalize_bin_transform_counts}
\title{Normalizes, bins and transforms counts}
\usage{
normalize_bin_transform_counts(
  count_matrix,
  idxstats_file,
  ref_size = 1e+08,
  bin_method = c("M5", "M24", "M12", "M3", "M2", "M1"),
  transform = c("asinh", "log2", "sqrt", "none")
)
}
\arguments{
\item{count_matrix}{DNase or ATAC-seq read counts matrix.}

\item{idxstats_file}{The \code{idxstats} file (generated by \code{samtools}).}

\item{ref_size}{Scale to DNase-seq or ATAC-seq reference library size.
(Default: 1e8 for DNase-seq and 5e7 for ATAC-seq).}

\item{bin_method}{\code{MILLIPEDE} binning scheme (Default: \sQuote{M5}).}

\item{transform}{Type of transformation for DNase or ATAC counts.
Options: \sQuote{asinh}, \sQuote{log2}, \sQuote{sqrt}, \sQuote{none}.}
}
\value{
A data frame of normalized, binned and transformed counts.
}
\description{
Normalizes counts by library size,
bin using \code{MILLIPEDE} binning method and then take \sQuote{asinh} or \sQuote{log2} transform.
}
\examples{
\dontrun{
# Normalizes counts by scaling to a library size of 100 million reads,
# and bin counts using MILLIPEDE M5 binning method, and then takes
# asinh transformation on the binned counts.
bins <- normalize_bin_transform_counts(count_matrix,
                                       idxstats_file,
                                       ref_size = 1e8,
                                       bin_method = 'M5',
                                       transform = 'asinh')
}
}
