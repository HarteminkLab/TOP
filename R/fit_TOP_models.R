
#' @title Fit TOP quantitative occupancy model with M5 bins using JAGS
#'
#' @param data a data frame containing the combined training data.
#' @param model.file file containing the TOP model written in BUGS code.
#' @param transform Type of transformation for ChIP counts.
#' Possible values are "asinh", "log2", "sqrt", and "none" (no transformation).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.chains number of Markov chains (default: 3).
#' @param n.thin thinning rate, must be a positive integer.
#' Default is max(1, floor(n.chains * (n.iter-n.burnin) / 1000))
#' which will only thin if there are at least 2000 simulations.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#' @export
fit_TOP_M5_model_jags <- function(data,
                                  model.file,
                                  transform=c('asinh', 'log2', 'sqrt', 'none'),
                                  n.iter=2000,
                                  n.burnin=floor(n.iter/2),
                                  n.chains=3,
                                  n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                                  quiet = FALSE) {

  if (!requireNamespace("R2jags", quietly = TRUE)) {
    stop(
      "Package \"R2jags\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(!all(c('pwm',paste0('bin', 1:5),'chip','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }

  transform <- match.arg(transform)

  # Transform the ChIP counts
  if (transform == 'asinh') {
    data$chip <- asinh(data$chip)
  } else if (transform == 'log2') {
    data$chip <- log2(data$chip + 1)
  } else if (transform == 'sqrt') {
    data$chip <- sqrt(data$chip)
  }else{
    cat('No transform done for ChIP data.\n')
  }

  # Training data
  training.data <- list(pwm = data$pwm,
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
  model.params <- c('A', 'Alpha', 'alpha',
                    'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
                    'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                    'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6',
                    'T', 'Tau', 'tau')

  # Fit Top M5 model
  jagsfit <- R2jags::jags(data = training.data,
                          parameters.to.save = model.params,
                          model.file = model.file,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          n.chains = n.chains,
                          quiet = quiet)

  return(jagsfit)

}


#' @title Fit TOP binary (logistic) model with M5 bins using JAGS
#'
#' @param data a data frame containing the combined training data.
#' @param model.file file containing the TOP model written in BUGS code.
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.chains number of Markov chains (default: 3).
#' @param n.thin thinning rate, must be a positive integer.
#' Default is max(1, floor(n.chains * (n.iter-n.burnin) / 1000))
#' which will only thin if there are at least 2000 simulations.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#' @export
#'
fit_TOP_logistic_M5_model_jags <- function(data,
                                           model.file,
                                           n.iter=2000,
                                           n.burnin=floor(n.iter/2),
                                           n.chains=3,
                                           n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                                           quiet = FALSE) {

  if (!requireNamespace("R2jags", quietly = TRUE)) {
    stop(
      "Package \"R2jags\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(!all(c('pwm',paste0('bin', 1:5),'chip_label','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }

  # Training data
  training.data <- list(pwm = data$pwm,
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
  model.params <- c('A', 'Alpha','alpha',
                    'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
                    'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                    'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6')

  ## Train Top M5 model using R2jags
  jagsfit <- R2jags::jags(data = training.data,
                          parameters.to.save = model.params,
                          model.file = model.file,
                          n.iter = n.iter,
                          n.burnin = n.burnin,
                          n.thin = n.thin,
                          n.chains = n.chains,
                          quiet = quiet)

  return(jagsfit)

}

#' @title a wrapper function to fit TOP model for the selected partitions in parallel
#'
#' @param all_training_data a list of the assembled training data of all partitions.
#' @param all_training_data_files a vector of the assembled training data
#' files of all partitions. If all_training_data is missing,
#' it will load the training data from all_training_data_files.
#' @param model.file file containing the TOP model written in BUGS code.
#' @param logistic.model Logical; if TRUE, use the logistic version of TOP model.
#' @param out.dir Output directory for TOP model posterior samples.
#' @param transform Type of transformation for ChIP counts.
#' Possible values are "asinh", "log2", "sqrt", and "none" (no transformation).
#' Only needed when logistic.model is FALSE.
#' @param partitions a vector selecting which partition(s) to run.
#' (default: all 10 partitions (1:10))
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.chains number of Markov chains (default: 3).
#' @param n.thin thinning rate, must be a positive integer.
#' Default is max(1, floor(n.chains * (n.iter-n.burnin) / 1000))
#' which will only thin if there are at least 2000 simulations.
#' @param n.cores Number of cores to use in parallel
#' (default: equal to the number of partitions).
#' @param quiet Logical, whether to suppress stdout in jags.model().
#'
#' @import doParallel
#' @import foreach
#'
#' @export
#'
fit_TOP_model <- function(all_training_data,
                          all_training_data_files,
                          model.file,
                          logistic.model = FALSE,
                          out.dir = "TOP_samples",
                          transform = c('asinh', 'log2', 'sqrt', 'none'),
                          partitions=1:10,
                          n.iter=2000,
                          n.burnin=floor(n.iter/2),
                          n.chains=3,
                          n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                          n.cores = length(partitions),
                          quiet = FALSE){

  if(!dir.exists(out.dir)){
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  }

  cat('Fitting TOP models for partition:', partitions,'...\n')

  registerDoParallel(cores=n.cores)
  cat('Using', getDoParWorkers(), 'cores. \n')

  # We can also submit jobs for each of the partitions
  # on separate nodes on compute clusters.
  TOP_samples_files <- foreach(k=partitions, .combine = "rbind") %dopar% {

    if(length(all_training_data) == 10){
      data <- all_training_data[[k]]
    }else if (length(all_training_data_files) == 10) {
      data <- readRDS(all_training_data_files[k])
    }else{
      stop('Check all_training_data or all_training_data_files!')
    }

    if(logistic.model){
      # obtain TOP logistic model posterior samples
      jagsfit <- fit_TOP_logistic_M5_model_jags(data, model.file, n.iter, n.burnin, n.chains, n.thin, quiet)
      fit_samples <- coda::as.mcmc(jagsfit)
      saveRDS(jagsfit, file.path(out.dir, paste0('TOP_logistic_M5_partition', k, '.jagsfit.rds')))
      samples_file <- file.path(out.dir, paste0('TOP_logistic_M5_partition', k, '.posterior_samples.rds'))
      saveRDS(fit_samples, samples_file)
    }else{
      # obtain TOP model posterior samples
      jagsfit <- fit_TOP_M5_model_jags(data, model.file, transform, n.iter, n.burnin, n.chains, n.thin, quiet)
      fit_samples <- coda::as.mcmc(jagsfit)
      saveRDS(jagsfit, file.path(out.dir, paste0('TOP_M5_partition', k, '.jagsfit.rds')))
      samples_file <- file.path(out.dir, paste0('TOP_M5_partition', k, '.posterior_samples.rds'))
      saveRDS(fit_samples, samples_file)
    }
    samples_file
  }

  cat("Files of TOP posterior samples: \n")
  print(TOP_samples_files)
}
