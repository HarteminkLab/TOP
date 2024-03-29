% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_chip_data.R
\name{normalize_chip}
\alias{normalize_chip}
\title{Normalizes (and transforms) ChIP-seq read coverage}
\usage{
normalize_chip(
  chip_counts,
  idxstats_file,
  ref_size = 2e+07,
  transform = c("none", "asinh", "log2", "sqrt")
)
}
\arguments{
\item{chip_counts}{A vector of ChIP-seq counts.}

\item{idxstats_file}{The \code{idxstats} file generated by \code{samtools}.}

\item{ref_size}{ChIP-Seq reference library size (Default: 20 million).}

\item{transform}{Type of transformation for the ChIP read counts.
Options are \sQuote{none} (no transformation), \sQuote{asinh}, \sQuote{log2}, and \sQuote{'sqrt'}.}
}
\value{
A vector of transformed ChIP-seq counts.
}
\description{
Normalizes the ChIP-seq read counts by scaling
to the reference ChIP-seq library size,
and transforms the ChIP-seq read counts.
}
