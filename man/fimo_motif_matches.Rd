% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_motif_sites.R
\name{fimo_motif_matches}
\alias{fimo_motif_matches}
\title{Runs \code{FIMO} to scan for motif matches}
\usage{
fimo_motif_matches(
  motif_file,
  sequence_file,
  outname = "fimo.txt",
  outdir = dirname(outname),
  thresh_pValue = 1e-05,
  background = c("default", "motif", "uniform", "file"),
  background_file,
  skip_matched_sequence = TRUE,
  max_strand = FALSE,
  max_stored_scores = 1e+05,
  options = "",
  verbosity = 2,
  fimo_path = "fimo"
)
}
\arguments{
\item{motif_file}{Motif file in MEME format.}

\item{sequence_file}{Genome sequence file in FASTA format.}

\item{outname}{Output file name.}

\item{outdir}{Output directory.}

\item{thresh_pValue}{\code{FIMO} option \code{thresh} for p-value threshold.}

\item{background}{Option for background model:
\sQuote{default}: uses \code{FIMO} default background setting;
\sQuote{motif}: uses the 0-order letter frequencies contained in the motif file;
\sQuote{uniform}: uses uniform letter frequencies;
\sQuote{file}: uses the file specified in \code{background_file}.}

\item{background_file}{Path to a file in Markov Background Model Format.}

\item{skip_matched_sequence}{\code{FIMO} option \code{skip_matched_sequence}.
Turns off output of the sequence of motif matches.
This speeds up processing considerably.}

\item{max_strand}{\code{FIMO} option \code{max_strand}.
If matches on both strands at a
given position satisfy the output threshold,
only report the match for the strand with the higher score.
If the scores are tied, the matching strand is chosen at random.}

\item{max_stored_scores}{The maximum number of stored matches.}

\item{options}{Other options for \code{FIMO}.}

\item{verbosity}{A number of the verbosity level (from 1 to 5).
If set to 1 (quiet) then it will only output error messages, in contrast, the
other extreme 5 (dump) outputs lots of mostly useless information.}

\item{fimo_path}{Path to \code{fimo} command line executable.
(default: 'fimo')}
}
\description{
Runs \code{FIMO} to scan for motif matches along the genome.
}
\examples{
\dontrun{
fimo_motif_matches(motif_file='motifID.meme',
                   sequence_file='hg38.fa',
                   thresh_pValue=1e-5,
                   outname='motifID_1e-5.fimo.txt',
                   fimo_path='fimo')
}
}
