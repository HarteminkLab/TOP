
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

    data <- read.table(data_file, header = TRUE, stringsAsFactors = FALSE)

    part_length <- as.integer(nrow(data) / total_partition) + 1
    selection <- (1 + (n_part - 1) * part_length):min(nrow(data), n_part*part_length)

    data <- data[selection, ]
    cat(tf_name, cell_type, ':', nrow(data), 'sites \n')
    data$tf <- tf_name
    data$cell_type <- cell_type
    combined_data[[ paste(tf_name, cell_type, sep = '|') ]] <- data
    n_data_read = n_data_read + 1

  }

  cat('Loaded', n_data_read, 'datasets.\n')

  ## row combine all data
  combined_data <- do.call(rbind, combined_data)
  row.names(combined_data) <- NULL

  return(combined_data)

}

#' @title Fit TOP model
#'
#' @param data combined training data.
#' @param model TOP model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#' @param DIC logical; if TRUE (default), compute deviance, pD, and DIC.
#'
#' @importFrom R2jags jags
#'
#' @export
#'
fit_TOP_M5_model <- function(data, model,
                             n.iter=10000, n.burnin=5000, n.thin=10, n.chains=1,
                             DIC=TRUE) {

  if(!all(c('pmw',paste0('bin', 1:5),'chip','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }
  # Check column names
  # colnames(data)[grep('pwm', colnames(data), ignore.case = T)] <- 'pwm'
  # colnames(data)[grep('dnase|atac', colnames(data), ignore.case = T)] <- paste0('bin', 1:5)
  # colnames(data)[grep('chip', colnames(data), ignore.case = T)] <- 'chip'

  # Training data
  training_data <- list(pwm = data$pwm,
                        bin1 = data$bin1,
                        bin2 = data$bin2,
                        bin3 = data$bin3,
                        bin4 = data$bin4,
                        bin5 = data$bin5,
                        chip = data$chip,
                        tf = data$tf_id,
                        cell_type = data$cell_id,
                        N = nrow(data),
                        n_tfs = length(unique(data$tf_id)),
                        n_cell_types = length(unique(data$cell_id)))

  # List parameters to be monitored
  model_params <- c('A', 'Alpha', 'alpha',
                    'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
                    'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                    'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6',
                    'T', 'Tau', 'tau')

  ## Fit Top M5 model using R2jags
  cat('Fit TOP M5 model... \n')

  jags_fit <- jags(data = training_data, model_params, model,
                  n.iter, n.burnin, n.thin, n.chains, DIC)

  return(jags_fit)

}


#' @title Fit TOP logistic model
#'
#' @param data combined training data.
#' @param model TOP logistic model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#' @param DIC logical; if TRUE (default), compute deviance, pD, and DIC.
#'
#' @importFrom R2jags jags
#'
#' @export
#'
fit_TOP_logistic_M5_model <- function(data, model,
                                      n.iter=10000, n.burnin=5000, n.thin=10, n.chains=1,
                                      DIC=TRUE) {

  if(!all(c('pmw',paste0('bin', 1:5),'chip','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }
  # Check column names
  # colnames(data)[grep('pwm', colnames(data), ignore.case = T)] <- 'pwm'
  # colnames(data)[grep('dnase|atac', colnames(data), ignore.case = T)] <- paste0('bin', 1:5)
  # colnames(data)[grep('chip', colnames(data), ignore.case = T)] <- 'chip'

  # Training data
  training_data <- list(pwm = data$pwm,
                        bin1 = data$bin1,
                        bin2 = data$bin2,
                        bin3 = data$bin3,
                        bin4 = data$bin4,
                        bin5 = data$bin5,
                        chip_label = data$chip_label,
                        tf = data$tf_id,
                        cell_type = data$cell_id,
                        N = nrow(data),
                        n_tfs = length(unique(data$tf_id)),
                        n_cell_types = length(unique(data$cell_id)))

  # List parameters to be monitored
  model_params <- c('A', 'Alpha','alpha',
                    'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
                    'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                    'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6')

  ## Fit Top M5 model using R2jags
  cat('Fit TOP logistic M5 model... \n')

  jags_fit <- jags(data = training_data, model_params, model,
                   n.iter, n.burnin, n.thin, n.chains, DIC)

  return(jags_fit)

}
