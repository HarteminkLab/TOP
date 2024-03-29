% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_profile}
\alias{plot_profile}
\title{Plots DNase or ATAC profiles}
\usage{
plot_profile(cuts, mlen = ncol(cuts)/2 - 200, title = "")
}
\arguments{
\item{cuts}{Matrix of DNase or ATAC cuts}

\item{mlen}{Motif length}

\item{title}{Title of the plot}
}
\description{
Plots DNase or ATAC profiles
}
\examples{
\dontrun{
plot_profile(cuts, title='Footprint profiles around binding sites')
}
}
