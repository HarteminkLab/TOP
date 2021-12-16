
#' @title Fit TOP model
#'
#' @param data combined training data.
#' @param model TOP model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#'
#' @importFrom R2jags jags
#'
#' @export
fit_TOP_M5_model <- function(data,
                             model.file,
                             n.iter=10000,
                             n.burnin=5000,
                             n.chains=3,
                             n.thin=10) {

  if(!all(c('pmw',paste0('bin', 1:5),'chip','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }

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

  jags_fit <- jags(data = training_data,
                   parameters.to.save = model_params,
                   model.file = model.file,
                   n.iter = n.iter,
                   n.burnin = n.burnin,
                   n.thin = n.thin,
                   n.chains = n.chains)

  return(jags_fit)

}


#' @title Fit TOP logistic model
#'
#' @param data combined training data.
#' @param model TOP logistic model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.chains number of Markov chains.
#' @param n.thin thinning rate, must be a positive integer.
#'
#' @importFrom R2jags jags
#'
#' @export
#'
fit_TOP_logistic_M5_model <- function(data,
                                      model.file,
                                      n.iter=10000,
                                      n.burnin=5000,
                                      n.chains=3,
                                      n.thin=10) {

  if(!all(c('pmw',paste0('bin', 1:5),'chip_label','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }

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

  jags_fit <- jags(data = training_data,
                   parameters.to.save = model_params,
                   model.file = model.file,
                   n.iter = n.iter,
                   n.burnin = n.burnin,
                   n.thin = n.thin,
                   n.chains = n.chains)

  return(jags_fit)

}
