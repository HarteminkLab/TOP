
#' @title Combine training data from TF x cell type combos
#'
#' @param tf_cell_table a data frame with TF name, PWM name, cell type, and link to the training data with PWM, DNase and ChIP data.
#' @param total_partition split the data into partitions (default = 10) and run Gibbs sampling in parallel.
#' @param n_part which partition to use for training the model
#'
#' @return The function returns a data frame of training data with all TFs and cell type combos.
#' @export
combine_training_data <- function(tf_cell_table, total_partition = 10, n_part = 1) {

  n_data_read = 0

  combined_data <- list()
  for(i in 1:nrow(tf_cell_table)) {

    tf_name <- as.character(tf_cell_table$tf_name[i])
    pwm_name <- as.character(tf_cell_table$pwm_name[i])
    cell_type <- as.character(tf_cell_table$cell_type[i])
    data_file <- as.character(tf_cell_table$data_file[i])

    if(!file.exists(data_file)){
      stop(data_file, 'not exist!')
    }

    data <- read.table(data_file, header = T, stringsAsFactors = F)

    part_length <- as.integer(nrow(data) / total_partition) + 1
    selection <- (1 + (n_part - 1) * part_length):min(nrow(data), n_part*part_length)

    data <- data[selection, ]
    cat(tf_name, cell_type, ':', nrow(data), 'sites \n')
    data$tf <- tf_name
    data$cell_type <- cell_type
    combined_data[[ paste(tf_name, cell_type, sep = '|') ]] <- data
    n_data_read = n_data_read + 1

  }

  cat('read', n_data_read, 'datasets\n')

  ## row combine all data
  combined_data <- do.call(rbind, combined_data)
  row.names(combined_data) <- NULL

  return(combined_data)

}

#' @title Fit TOP model using JAGS with R2jags package
#'
#' @param data.train combined training data.
#' @param TOP.model TOP model written in BUGS code.
#' @param parameters.to.save character vector of the names of the parameters to
#' save which should be monitored.
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#' @param DIC logical; if TRUE (default), compute deviance, pD, and DIC.
#'
#' @return
#' @export
fit_TOP_jags <- function(data.train, TOP.model, parameters.to.save,
                         n.iter=1e5, n.burnin=5000, n.thin=2, n.chains=1, DIC=TRUE) {

  attach(data.train)

  N <- nrow(data.train)
  n_tfs <- length(unique(data.train$tf_id))
  n_cell_types <- length(unique(data.train$cell_type_id))

  jags.data <- list('pwm',
                    'dnase.left2_sum', 'dnase.left1_sum', 'dnase.motif_sum', 'dnase.right1_sum', 'dnase.right2_sum',
                    'chip_count',
                    'cell_type_id', 'tf_id',
                    'N', 'n_tfs', 'n_cell_types')

  ## using R2jags
  jags_fit <- jags(data = jags.data, parameters.to.save, TOP.model,
                  n.iter, n.burnin, n.thin, n.chains, DIC)

  return(jags_fit)

}


