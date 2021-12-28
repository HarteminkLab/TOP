
#' @title Predict quantitative TF occupancy or TF binding probability
#' @description Predict quantitative TF occupancy or TF binding probability
#' using TOP model trained from ChIP-seq read counts or binary labels.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param TOP_coef A data frame or list containing the posterior samples
#' or posterior mean of TOP regression coefficients.
#' @param logistic.model Logical; if TRUE, use the logistic version of TOP model.
#' @param posterior.option Method to predict occupancy.
#' 'samples': uses posterior samples,
#' 'mean': uses posterior mean of regression coefficients.
#' @param transform Type of transformation for ChIP counts.
#' Possible values are "asinh", "log2", "sqrt", and "none" (no transformation).
#' Only needed when logistic.model is FALSE.
#'
#' @return Returns a vector of predicted TF occupancy (posterior mean).
#' @export
#'
predict_TOP <- function(data,
                        TOP_coef,
                        logistic.model = FALSE,
                        posterior.option = c('mean', 'samples'),
                        transform = c('asinh', 'log2', 'log', 'none')){

  posterior.option <- match.arg(posterior.option)

  if(logistic.model == FALSE) {
    transform <- match.arg(transform)
    if ( posterior.option == 'mean' ) {
      predictions <- predict_TOP_mean_coef(data, TOP_coef, transform = transform)
    } else if ( posterior.option == 'samples' ) {
      predictions <- predict_TOP_samples(data, TOP_coef, transform = transform)
    } else{
      stop('posterior.option needs to be samples or mean!')
    }
  }else if (logistic.model == TRUE){
    predictions <- predict_TOP_logistic_mean_coef(data, TOP_coef)
  }

  return(predictions)

}


#' @title Predict TF occupancy using posterior samples of regression coefficients
#' @description Predict TF occupancy using posterior samples of TOP regression
#' coefficients.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param coef_samples TOP posterior samples.
#' @param use.posterior.mean Logical; if TRUE, uses the posterior mean of
#' regression coefficients to make predictions.
#' @param sample.predictions Logical; if TRUE, sample from posterior predictions
#' and then take the mean of posterior prediction samples.
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Options: asinh, log2, sqrt, none.
#'
#' @return Returns a vector of predicted TF occupancy (posterior mean).
#' @export
predict_TOP_samples <- function(data,
                                coef_samples,
                                use.posterior.mean = FALSE,
                                sample.predictions = TRUE,
                                transform = c('asinh', 'log2', 'sqrt', 'none')){

  cat('Predicting TF occupancy using TOP posterior samples...\n')
  transform <- match.arg(transform)

  features <- select_features(data)
  data.matrix <- as.matrix(data.frame(intercept = 1, features, check.names = FALSE))

  alpha_samples <- coef_samples$alpha_samples
  beta_samples  <- coef_samples$beta_samples
  tau_samples   <- coef_samples$tau_samples
  coefficients  <- t(cbind(alpha_samples, beta_samples))

  if(use.posterior.mean){
    coefficients <- apply(coefficients, 1, mean)
    means <- data.matrix %*% coefficients
    predictions <- as.numeric(means)
  }else{
    means <- data.matrix %*% coefficients

    if(sample.predictions){
      sds <- 1 / sqrt(tau_samples)
      n_data <- nrow(data.matrix)
      n_samples <- ncol(coefficients)
      predictions <- sapply(1:n_samples, function(x) stats::rnorm(n_data, mean = means[, x], sd = sds[x]))
      predictions <- apply(predictions, 1, mean)
    }else{
      predictions <- as.numeric(apply(means, 1, mean))
    }

  }

  # transform back to the original scale
  if(transform == 'asinh'){
    predictions <- sinh(predictions)
  }else if (transform == 'log2'){
    # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'sqrt'){
    # sqrt(y)
    predictions <- predictions^2
  }else if (transform == 'none'){
    # no transform
    # cat('No transform done. \n')
  }

  return(predictions)

}


#' @title Predict TF occupancy using posterior mean of regression coefficients
#' @description Predict TF occupancy using posterior mean of TOP regression
#' coefficients.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param mean_coef A numeric vector. The posterior mean of trained regression
#' coefficients, including the intercept and coefficients for PWM score and
#' DNase (or ATAC) bins.
#' length(mean_coef) should be equal to 1+ncol(data).
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Options: asinh, log2, sqrt, none.
#'
#' @return Returns a vector of predicted TF occupancy (posterior mean).
#'
#' @export
#'
predict_TOP_mean_coef <- function(data,
                                  mean_coef,
                                  transform = c('asinh', 'log2', 'log', 'none')){

  cat('Predicting TF occupancy using TOP posterior mean coefficients...\n')
  transform <- match.arg(transform)

  features <- select_features(data)

  if((ncol(features)+1) != length(mean_coef)){
    stop('The number of coefficients not equal to
         the number of features + intercept! Check input data!')
  }

  coefficients <- as.matrix(mean_coef, ncol = 1)

  data.matrix <- as.matrix(data.frame(intercept = 1, features, check.names = FALSE))

  predictions <- as.numeric(data.matrix %*% coefficients)

  # transform back to the original scale
  if(transform == 'asinh'){
    predictions <- sinh(predictions)
  }else if (transform == 'log2'){
    # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'sqrt'){
    # sqrt(y)
    predictions <- predictions^2
  }else if (transform == 'none'){
    # no transform
    # cat('No transform done. \n')
  }

  return(predictions)

}

#' @title Predict TF binding probability by TOP logistic model
#' @description Predict TF binding probability using posterior mean of the regression
#' coefficients trained from TOP logistic model.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param mean_coef A numeric vector. The posterior mean of trained regression
#' coefficients, including the intercept and coefficients for PWM score and
#' DNase (or ATAC) bins.
#' length(mean_coef) should be equal to 1+ncol(data).
#'
#' @return Returns a vector of predicted TF binding probabilities.
#'
#' @export
#'
predict_TOP_logistic_mean_coef <- function(data, mean_coef){

  features <- select_features(data)

  cat('Predicting TF binding probability using TOP logistic model posterior mean coefficients...\n')

  if((ncol(features)+1) != length(mean_coef)){
    stop('The number of coefficients not equal to
         the number of features + intercept! Check input data!')
  }

  coefficients <- as.matrix(mean_coef, ncol = 1)

  data.matrix <- as.matrix(data.frame(intercept = 1, features, check.names = FALSE))

  mu <- as.numeric(data.matrix %*% coefficients)

  p <- 1 / (1 + exp(-mu)) # inverse logit

  return(p)

}

#' Select PWM and DNase (or ATAC) bin features from the input data frame
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param pwm.name name (prefix) of the PMW score column.
#' @param bin.name name (prefix) of the DNase or ATAC bin columns.
#'
select_features <- function(data, pwm.name = 'pwm', bin.name = 'bin'){
  data <- as.data.frame(data)
  pwm.col <- grep(pwm.name, colnames(data), ignore.case = TRUE, value = TRUE)
  bin.cols <- grep(bin.name, colnames(data), ignore.case = TRUE, value = TRUE)
  features.cols <- c(pwm.col, bin.cols)
  cat('Select features:', features.cols, '\n')
  features <- data[, features.cols]
  return(features)
}
