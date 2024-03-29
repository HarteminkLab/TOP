% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_dnase_atac_data.R
\name{merge_normalize_counts}
\alias{merge_normalize_counts}
\title{Merges DNase or ATAC-seq counts from multiple replicates,
then normalizes merged counts.}
\usage{
merge_normalize_counts(counts_files, idxstats_files, ref_size = 1e+08)
}
\arguments{
\item{counts_files}{DNase or ATAC-seq read counts matrix files.}

\item{idxstats_files}{The \code{idxstats} files generated by \code{samtools}.}

\item{ref_size}{Scales to reference library size.
(Default: 1e8 for DNase-seq and 5e7 for ATAC-seq).}
}
\value{
A matrix of merged and normalized counts.
}
\description{
Merges DNase or ATAC-seq read counts from multiple replicate samples,
and normalizes the read counts by scaling to the reference library size.
}
