
#' @title Train TOP model
#'
#' @param data combined training data.
#' @param model.file TOP model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.thin thinning rate, must be a positive integer.
#' @param n.chains number of Markov chains.
#' @param transform Transformation of ChIP counts (asinh, log2, sqrt, or none).
#' @param quiet Logical, whether to suppress stdout.
#'
#' @export
fit_TOP_M5_model_jags <- function(data,
                                  model.file,
                                  n.iter=2000,
                                  n.burnin=floor(n.iter/2),
                                  n.chains=3,
                                  n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                                  transform=c('asinh', 'log2', 'sqrt', 'none'),
                                  quiet = FALSE) {


  if(!all(c('pwm',paste0('bin', 1:5),'chip','tf_id','cell_id') %in% colnames(data))){
    stop('Check colnames of the data! \n')
  }

  transform <- match.arg(transform)

  # Transform the ChIP counts
  if (transform == 'asinh') {
    cat(transform, 'transform ChIP data...\n')
    data$chip <- asinh(data$chip)
  } else if (transform == 'log2') {
    cat(transform, 'transform ChIP data...\n')
    data$chip <- log2(data$chip + 1)
  } else if (transform == 'sqrt') {
    cat(transform, 'transform ChIP data...\n')
    data$chip <- sqrt(data$chip)
  }else{
    cat('No transform for ChIP data...\n')
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

  ## Fit Top M5 model using R2jags
  cat('Fit TOP M5 model... \n')

  if (requireNamespace("R2jags", quietly = TRUE)) {
    # Use "R2jags" package
    cat('using R2jags ...\n')
    jagsfit <- R2jags::jags(data = training.data,
                            parameters.to.save = model.params,
                            model.file = model.file,
                            n.iter = n.iter,
                            n.burnin = n.burnin,
                            n.thin = n.thin,
                            n.chains = n.chains,
                            quiet = quiet)
    fit.samples <- as.mcmc(jagsfit)
    cat('Done. \n')
  } else if (requireNamespace("rjags", quietly = TRUE)) {
    # Use "rjags" package
    # Set up the JAGS model and settings
    cat('using rjags ...\n')
    cat('Init model... \n')
    jags.obj <- rjags::jags.model(data = training.data, file = model.file, n.chains = n.chains, quiet = quiet)
    # Burn-ins
    cat(n.burnin, 'burn-in... \n')
    stats::update(jags.obj, n.iter = n.burnin)
    # run JAGS and save posterior samples as a mcmc.list object
    cat('Run Gibbs sampling... \n')
    fit.samples <- rjags::coda.samples(jags.obj, variable.names = model.params, n.iter = n.iter-n.burnin, thin = n.thin)
    cat('Done. \n')
  }else{
    stop(
      "Package \"R2jags\" or \"rjags\" must be installed to fit TOP model.",
      call. = FALSE
    )
  }

  return(fit.samples)

}


#' @title Train TOP logistic model
#'
#' @param data combined training data.
#' @param model.file TOP logistic model (written in BUGS code).
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.chains number of Markov chains.
#' @param n.thin thinning rate, must be a positive integer.
#' @param quiet Logical, whether to suppress stdout.
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
  cat('Train TOP logistic M5 model... \n')

  if (requireNamespace("R2jags", quietly = TRUE)) {
    # Use "R2jags" package
    cat('using R2jags ...\n')
    jagsfit <- R2jags::jags(data = training.data,
                            parameters.to.save = model.params,
                            model.file = model.file,
                            n.iter = n.iter,
                            n.burnin = n.burnin,
                            n.thin = n.thin,
                            n.chains = n.chains,
                            quiet = quiet)
    fit.samples <- as.mcmc(jagsfit)
    cat('Done. \n')
  } else if (requireNamespace("rjags", quietly = TRUE)) {
    # Use "rjags" package
    # Set up the JAGS model and settings
    cat('using rjags ...\n')
    cat('Init model... \n')
    jags.obj <- rjags::jags.model(data = training.data, file = model.file, n.chains = n.chains, quiet = quiet)
    # Burn-ins
    cat(n.burnin, 'burn-in... \n')
    stats::update(jags.obj, n.iter = n.burnin)
    # run JAGS and save posterior samples as a mcmc.list object
    cat('Run Gibbs sampling... \n')
    fit.samples <- rjags::coda.samples(jags.obj, variable.names = model.params, n.iter = n.iter-n.burnin, thin = n.thin)
    cat('Done. \n')
  }else{
    stop(
      "Package \"R2jags\" or \"rjags\" must be installed to fit TOP model.",
      call. = FALSE
    )
  }

  return(fit.samples)

}

#' Fit TOP model for each partition separately
#' @param data.partitions Training data of all partitions
#' @param model.file TOP logistic model file.
#' @param logistic.model If TRUE, use logistic version of the model
#' @param out.dir Output directory for TOP model posterior samples
#' @param partitions select which partition(s) to run
#' @param n.iter number of total iterations per chain (including burn in).
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning.
#' @param n.chains number of Markov chains.
#' @param n.thin thinning rate, must be a positive integer.
#' @param transform Transformation of ChIP counts (asinh, log2, sqrt, or none)
#' @param quiet Logical, whether to suppress stdout.
#' @import doParallel
#' @import foreach
#'
#' @export
#'
fit_TOP_model <- function(all.training.data,
                          all.training.data.files,
                          model.file,
                          logistic.model = FALSE,
                          out.dir = "TOP_samples",
                          partitions=1:10,
                          n.iter=2000,
                          n.burnin=floor(n.iter/2),
                          n.chains=3,
                          n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                          transform = c('asinh', 'log2', 'sqrt', 'none'),
                          quiet = TRUE,
                          parallel = TRUE){

  transform <- match.arg(transform)

  if(!dir.exists(out.dir)){
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  }

  partitions <- as.integer(partitions)

  if(parallel){
    doParallel::registerDoParallel(cores=length(partitions))
    cat("Using", foreach::getDoParWorkers(), "cores in parallel. \n")
  }

  # We can submit jobs in parallel for the partitions on separate compute nodes
  # instead of using the foreach loop here.
  TOP_samples_files <- foreach(k=partitions, .combine = "rbind") %dopar% {
    cat('Load assembled training data in partition: ', k, '\n')

    if(length(all.training.data) == 10){
      data <- all.training.data[[k]]
    }else if (length(all.training.data.files) == 10) {
      data <- readRDS(all.training.data.files[k])
    }else{
      stop('Check all.training.data or all.training.data.files!')
    }

    cat('Training TFs: ', unique(data$tf_name), '\n')
    cat('Training cell types: ', unique(data$cell_type), '\n')

    if(logistic.model){
      # obtain TOP logistic model posterior samples
      TOP_samples <- fit_TOP_logistic_M5_model_jags(data, model.file, n.iter, n.burnin, n.chains, n.thin, quiet)
      out_file <- paste0(out.dir, '/TOP_logistic_M5_partition', k, '.posterior_samples.rds')
      saveRDS(TOP_samples, out_file)
    }else{
      # obtain TOP model posterior samples
      TOP_samples <- fit_TOP_M5_model_jags(data, model.file, n.iter, n.burnin, n.chains, n.thin, transform, quiet)
      out_file <- paste0(out.dir, '/TOP_M5_partition', k, '.posterior_samples.rds')
      saveRDS(TOP_samples, out_file)
    }
    out_file
  }

  cat("Files of TOP model posterior samples: \n")
  print(TOP_samples_files)
}
