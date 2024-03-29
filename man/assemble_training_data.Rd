% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assemble_training_data.R
\name{assemble_training_data}
\alias{assemble_training_data}
\title{Assembles TOP training data for all TF x cell type combinations,
then split training data into 10 partitions}
\usage{
assemble_training_data(
  tf_cell_table,
  logistic_model = FALSE,
  chip_col = "chip",
  training_chrs = paste0("chr", seq(1, 21, 2)),
  n_partitions = 10,
  n_cores = n_partitions,
  max_sites = 50000,
  seed = 1
)
}
\arguments{
\item{tf_cell_table}{A data frame listing all TF x cell type combinations
and the training data for each combination.
It should have at least three columns, with:
TF names, cell types, and file names of the individual training data for
each TF x cell type combination. The individual training data
should be in .rds or text (.txt, or .csv) format.}

\item{logistic_model}{Logical. If \code{logistic_model = TRUE},
prepare assembled data for the logistic version of TOP model.
If \code{logistic_model = FALSE}, prepare assembled data for the
quantitative occupancy model (default).}

\item{chip_col}{The column name of ChIP data in the individual training data
(default: \sQuote{chip}).}

\item{training_chrs}{Chromosomes used for training the model
(default: odd chromosomes, chr1, chr3, ..., chr21)}

\item{n_partitions}{Number of partitions to split the training data (default: 10).}

\item{n_cores}{Number of cores to run in parallel
(default: equal to \code{n_partitions}).}

\item{max_sites}{Max number of candidate sites to keep for
each TF x cell type combination (default: 50000). To reduce computation time,
randomly select \code{max_sites} candidate sites for
each TF x cell type combination, if the number of candidate sites
exceeds \code{max_sites}.}

\item{seed}{A number for the seed used when sampling sites.}
}
\value{
A list of data frames (default: 10),
each containing one partition of the
training data with all TF x cell type combinations.
}
\description{
Prepares the training data for fitting TOP models.
It splits training data into 10 partitions and
assembles training data for all TF x cell type combinations
for each of the partitions.
}
\examples{
\dontrun{

#  tf_cell_table should have three columns with:
#  TF names, cell types, and paths to the training data files, like:
#  |   tf_name    |   cell_type   |        data_file         |
#  |:------------:|:-------------:|:------------------------:|
#  |     CTCF     |     K562      |   CTCF.K562.data.rds     |
#  |     CTCF     |     A549      |   CTCF.A549.data.rds     |
#  |     CTCF     |    GM12878    |   CTCF.GM12878.data.rds  |
#  |     ...      |     ...       |   ...                    |

# Assembles training data for the quantitative occupancy model,
# uses odd chromosomes for training, keeps at most 50000 candidate sites for
# each TF x cell type combination, and splits training data into 10 partitions.
assembled_training_data <- assemble_training_data(tf_cell_table,
                                                  logistic_model = FALSE,
                                                  chip_col = 'chip',
                                                  training_chrs = paste0('chr', seq(1,21,2)),
                                                  n_partitions=10,
                                                  max_sites = 50000)

# Assembles training data for the logistic version of the model
assembled_training_data <- assemble_training_data(tf_cell_table,
                                                  logistic_model = TRUE,
                                                  chip_col = 'chip_label',
                                                  training_chrs = paste0('chr', seq(1,21,2)),
                                                  n_partitions=10,
                                                  max_sites = 50000)

}

}
