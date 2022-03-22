
#' @title Fit TOP model
#' @description Fit TOP model (using M5 bins) for selected partitions in parallel.
#'
#' @param all_training_data A list of the assembled training data of all partitions.
#' @param all_training_data_files A vector of the assembled training data
#' files of all partitions. If all_training_data is missing,
#' it will load the training data from all_training_data_files.
#' @param model_file TOP model file written in JAGS.
#' @param logistic_model Logical; if TRUE, use the logistic version of TOP model.
#' @param transform Type of transformation for ChIP counts.
#' Possible values are "asinh", "log2", "sqrt", and "none" (no transformation).
#' Only needed when logistic_model is FALSE.
#' @param partitions A vector of selected partition(s) to run.
#' (default: all 10 partitions (1:10))
#' @param n_iter Number of total iterations per chain (including burn in).
#' @param n_burnin Length of burn in samples,
#' i.e. number of iterations to discard at the beginning.
#' Default is n_iter/2, that is, discarding the first half of the simulations.
#' @param n_chains Number of Markov chains (default: 3).
#' @param n_thin Thinning rate, must be a positive integer.
#' Default is max(1, floor(n_chains * (n_iter-n_burnin) / 1000))
#' which will only thin if there are at least 2000 simulations.
#' @param n_cores Number of cores to use in parallel
#' (default: equal to the number of partitions).
#' @param save Logical, if TRUE, save posterior samples as .rds files
#' in \code{out_dir}.
#' @param out_dir Directory to save TOP model posterior samples.
#' @param return_type Type of result to return.
#' Options: 'samples' (posterior samples),
# 'jagsfit' (jagsfit object), or 'none' (no return values).
#' @param quiet Logical, whether to suppress stdout in jags.model().
#' @return A list of posterior samples for each of the partitions.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores
#' @return A list of posterior samples or jagsfit object for each partition.
#' @export
#'
fit_TOP_M5_model <- function(all_training_data,
                             all_training_data_files,
                             model_file,
                             logistic_model=FALSE,
                             transform=c('asinh', 'log2', 'sqrt', 'none'),
                             partitions=1:10,
                             n_iter=2000,
                             n_burnin=floor(n_iter/2),
                             n_chains=3,
                             n_thin=max(1, floor((n_iter - n_burnin) / 1000)),
                             n_cores=length(partitions),
                             save=TRUE,
                             out_dir='./TOP_samples',
                             return_type=c('samples', 'jagsfit', 'samplefiles'),
                             quiet=FALSE){

  transform <- match.arg(transform)
  return_type <- match.arg(return_type)

  cat('Fitting TOP models for partition:', partitions,'...\n')

  n_partitions <- length(partitions)
  if(missing(n_cores)){
    n_available_cores <- detectCores(logical = FALSE) - 1
    n_cores <- min(n_available_cores, n_partitions)
  }else{
    n_cores <- min(n_cores, n_partitions)
  }

  registerDoParallel(cores=n_cores)
  cat('Using', getDoParWorkers(), 'cores. \n')

  # We can also run each of the partitions
  # on separate compute nodes.
  TOP_samples <- foreach(k=partitions, .combine = "rbind") %dopar% {

    if(length(all_training_data) == 10){
      data <- all_training_data[[k]]
    }else if (length(all_training_data_files) == 10) {
      data <- readRDS(all_training_data_files[k])
    }else{
      stop('Check training data!')
    }

    if(logistic_model){
      # Fit TOP binding probability model (logistic version)
      jagsfit <- fit_TOP_logistic_M5_model_jags(data, model_file, n_iter, n_burnin, n_chains, n_thin, quiet)
      if(save){
        if(!dir.exists(out_dir)) dir.create(out_dir)
        saveRDS(jagsfit, file.path(out_dir, paste0('TOP.logistic.M5.partition', k, '.jagsfit.rds')))
        samples_file <- file.path(out_dir, paste0('TOP.logistic.M5.partition', k, '.posterior.samples.rds'))
        saveRDS(coda::as.mcmc(jagsfit), samples_file)
      }

    }else{
      # Fit TOP quantitative TF occupancy model
      jagsfit <- fit_TOP_M5_model_jags(data, model_file, transform, n_iter, n_burnin, n_chains, n_thin, quiet)
      if(save){
        if(!dir.exists(out_dir)) dir.create(out_dir)
        saveRDS(jagsfit, file.path(out_dir, paste0('TOP.M5.partition', k, '.jagsfit.rds')))
        samples_file <- file.path(out_dir, paste0('TOP.M5.partition', k, '.posterior.samples.rds'))
        saveRDS(coda::as.mcmc(jagsfit), samples_file)
      }
    }

    if(return_type == 'samples'){
      coda::as.mcmc(jagsfit)
    }else if(return_type == 'jagsfit'){
      jagsfit
    }else if(return_type == 'samplefiles' && save){
      samples_file
    }

  }

  return(TOP_samples)
}

#' @title Fit TOP quantitative occupancy model with M5 bins
#' @description Fit TOP quantitative occupancy model with M5 bins using JAGS.
#' The R2jags package \code{\link[R2jags]{R2jags}} is required.

#' @param data a data frame containing the combined training data.
#' @param model_file TOP model file written in JAGS.
#' @param transform Type of transformation for ChIP counts.
#' Possible values are "asinh", "log2", "sqrt", and "none" (no transformation).
#' @param n_iter number of total iterations per chain (including burn in).
#' @param n_burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' Default is n_iter/2, that is, discarding the first half of the simulations.
#' @param n_chains number of Markov chains (default: 3).
#' @param n_thin thinning rate, must be a positive integer.
#' Default is max(1, floor(n_chains * (n_iter-n_burnin) / 1000))
#' which will only thin if there are at least 2000 simulations.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#' @return A jagsfit object from \code{R2jags}.
#' @export
fit_TOP_M5_model_jags <- function(data,
                                  model_file,
                                  transform=c('asinh', 'log2', 'sqrt', 'none'),
                                  n_iter=2000,
                                  n_burnin=floor(n_iter/2),
                                  n_chains=3,
                                  n_thin=max(1, floor((n_iter - n_burnin) / 1000)),
                                  quiet=FALSE) {

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

  # Transform ChIP counts
  if (transform == 'asinh') {
    data$chip <- asinh(data$chip)
  } else if (transform == 'log2') {
    data$chip <- log2(data$chip + 1)
  } else if (transform == 'sqrt') {
    data$chip <- sqrt(data$chip)
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

  # List parameters to save
  model_params <- c('A', 'Alpha', 'alpha',
                    'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
                    'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                    'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6',
                    'T', 'Tau', 'tau')

  # Fit Top M5 model using R2jags
  jagsfit <- R2jags::jags(data = training_data,
                          parameters.to.save = model_params,
                          model.file = model_file,
                          n.iter = n_iter,
                          n.burnin = n_burnin,
                          n.thin = n_thin,
                          n.chains = n_chains,
                          quiet = quiet)

  return(jagsfit)

}


#' @title Fit TOP logistic model with M5 DNase (or ATAC) bins and binary ChIP labels.
#' @description Fit TOP logistic model with M5 DNase (or ATAC) bins and
#' binary ChIP labels using JAGS.
#' The R2jags package \code{\link[R2jags]{R2jags}} is required.
#' @param data a data frame containing the combined training data.
#' @param model_file TOP model file written in JAGS.
#' @param n_iter number of total iterations per chain (including burn in).
#' @param n_burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' Default is n_iter/2, that is, discarding the first half of the simulations.
#' @param n_chains number of Markov chains (default: 3).
#' @param n_thin thinning rate, must be a positive integer.
#' Default is max(1, floor(n_chains * (n_iter-n_burnin) / 1000))
#' which will only thin if there are at least 2000 simulations.
#' @param quiet Logical, whether to suppress stdout in jags.model().
#' @return A jagsfit object from \code{R2jags}.
#' @export
#'
fit_TOP_logistic_M5_model_jags <- function(data,
                                           model_file,
                                           n_iter=2000,
                                           n_burnin=floor(n_iter/2),
                                           n_chains=3,
                                           n_thin=max(1, floor((n_iter - n_burnin) / 1000)),
                                           quiet=FALSE) {

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

  # List parameters to save
  model_params <- c('A', 'Alpha','alpha',
                    'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
                    'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
                    'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6')

  # Fit Top M5 logistic model using R2jags
  jagsfit <- R2jags::jags(data = training_data,
                          parameters.to.save = model_params,
                          model.file = model_file,
                          n.iter = n_iter,
                          n.burnin = n_burnin,
                          n.thin = n_thin,
                          n.chains = n_chains,
                          quiet = quiet)

  return(jagsfit)

}
