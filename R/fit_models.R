
#' @title Fit TOP model with M5 bins
#' @description Fit TOP model with M5 bins for selected partitions in parallel.
#' By default, it runs Gibbs sampling for all 10 partitions
#' in parallel on 10 CPU cores, and returns a list of posterior samples
#' for each of the 10 partitions.
#' Alternatively, you may fit model for each of the 10
#' the partitions on separate machines by specifying the partition to run.
#' @param all_training_data A list of the assembled training data of all partitions.
#' @param all_training_data_files A vector of the assembled training data
#' files of all partitions. If \code{all_training_data} is missing,
#' it will load the training data from \code{all_training_data_files}.
#' @param model_file TOP model file written in \code{JAGS}.
#' @param logistic_model Logical; whether to use the logistic version of TOP model.
#' If \code{logistic_model = TRUE}, use the logistic version of TOP model.
#' If \code{logistic_model = FALSE}, use the quantitative occupancy model (default).
#' @param transform Type of transformation for ChIP-seq read counts.
#' Options are: \sQuote{asinh}(asinh transformation),
#' \sQuote{log2} (log2 transformation),
#' \sQuote{sqrt} (square root transformation),
#' and \sQuote{none}(no transformation).
#' This only applies when \code{logistic_model = FALSE}.
#' @param partitions A vector of selected partition(s) to run.
#' Default: all 10 partitions. If you specify a few partitions,
#' it will only fit models to data in those selected partitions.
#' @param n_iter Number of total iterations per chain, including burn-in iterations.
#' @param n_burnin Length of burn-in iterations,
#' i.e. number of samples to discard at the beginning.
#' Default is \code{n_iter/2}, discarding the first half of the samples.
#' @param n_chains Number of Markov chains (default: 3).
#' @param n_thin Thinning rate, must be a positive integer.
#' Default is \code{max(1, floor(n_chains * (n_iter-n_burnin) / 1000))}
#' which will only thin if there are at least 2000 simulations.
#' No thinning will be performed if \code{n_thin = 1}.
#' @param n_cores Number of cores to use in parallel
#' (default: equal to the number of partitions, i.e. \code{length(partitions)}).
#' @param save Logical, if TRUE, save posterior samples as \sQuote{.rds} files
#' in \code{outdir}.
#' @param outdir Directory to save TOP model posterior samples.
#' @param return_type Specify the type of result to return.
#' Options: \sQuote{samples}(posterior samples),
#' \sQuote{jagsfit} (\code{jagsfit} object),
#' or \sQuote{samplefiles} (file names of posterior samples).
#' @param quiet Logical, if TRUE, suppress model fitting messages.
#' Otherwise, only show progress bars.
#' @return A list of posterior samples for each of the partitions.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores
#' @return A list of posterior samples or \code{jagsfit} object for each partition.
#' @export
#' @examples
#'
#' We can train the TOP model after we got the training data assembled.
#' Please read the data preparation tutorial to
#' prepare the training data.
#'
#' \dontrun{
#' # Example to train TOP quantitative occupancy model:
#' model_file <- system.file("model", "TOP_M5_model.jags", package = "TOP")
#'
#' # The example below first performs "asinh" transform to the ChIP-seq counts
#' # in 'assembled_training_data', then runs Gibbs sampling
#' # for each of the 10 partitions in parallel.
#' # The following example runs 5000 iterations of Gibbs sampling in total,
#' # including 1000 burn-ins, with 3 Markov chains, at a thinning rate of 2,
#' # and save the posterior samples to the TOP_fit directory.
#'
#' all_TOP_samples <- fit_TOP_M5_model(assembled_training_data,
#'                                     model_file = model_file,
#'                                     logistic_model = FALSE,
#'                                     transform = 'asinh',
#'                                     partitions = 1:10,
#'                                     n_iter = 5000,
#'                                     n_burnin = 1000,
#'                                     n_chains = 3,
#'                                     n_thin = 2,
#'                                     save = TRUE,
#'                                     out_dir = 'TOP_fit',
#'                                     quiet = TRUE)
#'
#' # We can also obtain the posterior samples separately for each partition,
#' # For example, to obtain the posterior samples for partition 3 only:
#'
#' TOP_samples_part3 <- fit_TOP_M5_model(assembled_training_data,
#'                                       model_file = model_file,
#'                                       logistic_model = FALSE,
#'                                       transform = 'asinh',
#'                                       partitions = 3,
#'                                       n_iter = 5000,
#'                                       n_burnin = 1000,
#'                                       n_chains = 3,
#'                                       n_thin = 2,
#'                                       save = TRUE,
#'                                       out_dir = 'TOP_fit',
#'                                       quiet = TRUE)
#'
#'
#' # Example to train TOP logistic (binary) model:
#' model_file <- system.file("model", "TOP_M5_logistic_model.jags", package = "TOP")
#'
#' all_TOP_samples <- fit_TOP_M5_model(assembled_training_data,
#'                                     model_file = model_file,
#'                                     logistic_model = TRUE,
#'                                     partitions = 1:10,
#'                                     n_iter = 5000,
#'                                     n_burnin = 1000,
#'                                     n_chains = 3,
#'                                     n_thin = 2,
#'                                     save = TRUE,
#'                                     out_dir = 'TOP_fit',
#'                                     quiet = TRUE)
#'
#' }
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
                             outdir='./TOP_samples',
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
      stop('Training data must contain 10 partitions!')
    }

    if(logistic_model){
      # Fit TOP binding probability model (logistic version)
      jagsfit <- fit_TOP_logistic_M5_model_jags(data, model_file, n_iter, n_burnin, n_chains, n_thin, quiet)
      if(save){
        if(!dir.exists(outdir)) dir.create(outdir)
        saveRDS(jagsfit, file.path(outdir, paste0('TOP.logistic.M5.partition', k, '.jagsfit.rds')))
        samples_file <- file.path(outdir, paste0('TOP.logistic.M5.partition', k, '.posterior.samples.rds'))
        saveRDS(coda::as.mcmc(jagsfit), samples_file)
      }

    }else{
      # Fit TOP quantitative occupancy model
      jagsfit <- fit_TOP_occupancy_M5_model_jags(data, model_file, transform, n_iter, n_burnin, n_chains, n_thin, quiet)
      if(save){
        if(!dir.exists(outdir)) dir.create(outdir)
        saveRDS(jagsfit, file.path(outdir, paste0('TOP.M5.partition', k, '.jagsfit.rds')))
        samples_file <- file.path(outdir, paste0('TOP.M5.partition', k, '.posterior.samples.rds'))
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

# Fit TOP quantitative occupancy model with M5 bins and quantitative ChIP occupancy (read counts)
# The R2jags package \code{\link[R2jags]{R2jags}} is required.
fit_TOP_occupancy_M5_model_jags <- function(data,
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

  if(anyNA(data)){
    stop('Data contains missing values! \n')
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


# Fit TOP logistic model with M5 DNase (or ATAC) bins and binary ChIP labels.
# The R2jags package \code{\link[R2jags]{R2jags}} is required.
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

  if(anyNA(data)){
    stop('Data contains missing values! \n')
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
