
#' @title Train TOP model
#'
#' @param data combined training data.
#' @param model TOP model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#'
#' @import R2jags
#'
#' @export
train_TOP_M5_model_jags <- function(data,
                               model,
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

  ## Train Top M5 model using R2jags
  cat('Train TOP M5 model... \n')

  TOP_samples <- R2jags::jags(data = training_data,
                              parameters.to.save = model_params,
                              model.file = model,
                              n.iter = n.iter,
                              n.burnin = n.burnin,
                              n.thin = n.thin,
                              n.chains = n.chains)

  return(TOP_samples)

}


#' @title Train TOP logistic model
#'
#' @param data combined training data.
#' @param model TOP logistic model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.chains number of Markov chains.
#' @param n.thin thinning rate, must be a positive integer.
#'
#' @import R2jags
#'
#' @export
#'
train_TOP_logistic_M5_model_jags <- function(data,
                                        model,
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

  ## Train Top M5 model using R2jags
  cat('Train TOP logistic M5 model... \n')

  TOP_samples <- R2jags::jags(data = training_data,
                              parameters.to.save = model_params,
                              model.file = model,
                              n.iter = n.iter,
                              n.burnin = n.burnin,
                              n.thin = n.thin,
                              n.chains = n.chains)

  return(TOP_samples)

}


#' Train TOP model for each partition separately
#' @param model_file TOP logistic model file.
#' @param training_data_dir Directory for saving training data
#' @param training_data_name Prefix for training data file names
#' @param logistic.model If TRUE, use logistic version of the model
#' @param out_dir Output directory for TOP model posterior samples
#' @param partitions select which partition(s) to run
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.chains number of Markov chains.
#' @param n.thin thinning rate, must be a positive integer.
#' @import doParallel
#' @import foreach
#'
#' @export
#'
train_TOP_model <- function(model_file,
                            training_data_dir,
                            training_data_name,
                            logistic.model = FALSE,
                            out_dir,
                            partitions=1:10,
                            n.iter=10000,
                            n.burnin=5000,
                            n.chains=3,
                            n.thin=10){

  if(!dir.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }

  partitions <- as.integer(partitions)

  doParallel::registerDoParallel(cores=length(partitions))
  cat("Using", foreach::getDoParWorkers(), "cores in parallel. \n")

  # We can submit jobs in parallel for the partitions on separate compute nodes
  # instead of using the foreach loop here.
  res_file_list <- foreach(k=partitions, .combine = "rbind") %dopar% {
    training_data <- readRDS(file.path(training_data_dir, paste0(training_data_name, '.partition', k, '.rds')))
    cat('Partition: ', k, '\n')
    cat('Training TFs: ', levels(training_data$tf_name), '\n')
    cat('Training cell types: ', levels(training_data$cell_type), '\n')

    if(logistic.model){
      # get TOP logistic model posterior samples
      TOP_samples <- train_TOP_logistic_M5_model_jags(training_data, model_file, n.iter, n.burnin, n.chains, n.thin)
      out_file <- paste0(out_dir, '/TOP_logistic_M5_partition', k, '.posterior_samples.rds')
      saveRDS(TOP_samples, out_file)
    }else{
      # get TOP model posterior samples
      TOP_samples <- train_TOP_M5_model_jags(training_data, model_file, n.iter, n.burnin, n.chains, n.thin)
      out_file <- paste0(out_dir, '/TOP_M5_partition', k, '.posterior_samples.rds')
      saveRDS(TOP_samples, out_file)
    }
    out_file
  }

  cat("Files of TOP model posterior samples: \n")
  print(res_file_list)
}
