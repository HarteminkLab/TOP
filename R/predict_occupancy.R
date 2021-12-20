
#' @title Predict quantitative TF occupancy or TF binding probability
#' @description Predict quantitative TF occupancy or TF binding probability
#' using TOP model trained from ChIP-seq read counts or binary labels.
#'
#' @param data A data frame. Columns are motif score and DNase (or ATAC) bins.
#' Rows are candidate sites.
#' @param TOP_model TOP posterior samples or posterior mean of regression
#' coefficients.
#' @param logistic.model If TRUE, use logistic version of the model
#' to predict TF binding probability.
#' @param posterior.option "samples": uses posterior samples, "mean": uses posterior mean of
#' trained regression coefficients
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP quantitative model. Options: asinh, log2, log, none.
#' This does not apply to the logistic version.
#'
#' @return The function returns a vector of predicted TF occupancy.
#' @export
predict_TOP <- function(data,
                        TOP_model,
                        logistic.model = FALSE,
                        posterior.option = c("mean", "samples"),
                        transform = c('asinh', 'log2', 'log', 'none')){

  posterior.option <- match.arg(posterior.option)

  if(logistic.model == FALSE) {
    transform <- match.arg(transform)
    if ( posterior.option == 'mean' ) {
      predictions <- predict_TOP_mean_coef(data, TOP_model, transform = transform)
    } else if ( posterior.option == 'samples' ) {
      predictions <- predict_TOP_samples(data, TOP_model, transform = transform)
    } else{
      stop('posterior.option needs to be samples or mean!')
    }
  }else if (logistic.model == TRUE){
    predictions <- predict_TOP_logistic_mean_coef(data, TOP_model)
  }

  return(predictions)

}


#' @title Predict TF occupancy using posterior samples of regression coefficients
#' @description Predict TF occupancy using posterior samples of regression
#' coefficients trained from the TOP model
#'
#' @param data A data frame. Columns are motif score and DNase (or ATAC) bins.
#' Rows are candidate sites.
#' @param TOP_samples TOP posterior samples.
#' @param use.posterior.mean If TRUE, uses the posterior mean of
#' regression coefficients to make predictions.
#' @param sample.predictions If TRUE, sample from posterior predictions and then
#' take the mean of posterior prediction samples.
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Options: asinh, log2, log, none.
#'
#' @return The function returns a vector of predicted TF occupancy.
#' @export
predict_TOP_samples <- function(data,
                                TOP_samples,
                                use.posterior.mean = FALSE,
                                sample.predictions = TRUE,
                                transform = c('asinh', 'log2', 'log', 'none')){

  transform <- match.arg(transform)

  alpha_samples <- TOP_samples$alpha_samples
  beta_samples  <- TOP_samples$beta_samples
  tau_samples   <- TOP_samples$tau_samples

  coefficients  <- t(cbind(alpha_samples, beta_samples))

  data <- as.matrix(data.frame(intercept = 1, data, check.names = F))

  if(use.posterior.mean){
    coefficients <- apply(coefficients, 1, mean)
    means <- data %*% coefficients
    predictions <- as.numeric(means)
  }else{
    means <- data %*% coefficients

    if(sample.predictions){
      sds <- 1 / sqrt(tau_samples)

      n_data <- nrow(data)
      n_samples <- ncol(coefficients)

      predictions <- sapply(1:n_samples, function(x) rnorm(n_data, mean = means[, x], sd = sds[x]))
      predictions <- apply(predictions, 1, mean)
    }else{
      predictions <- as.numeric(apply(means, 1, mean))
    }

  }

  if(transform == 'asinh'){
    predictions <- sinh(predictions)
  }else if (transform == 'log2'){
    # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'log'){
    # log(y+1)
    predictions <- exp(predictions) - 1
  }else if (transform == 'none'){
    # no transform
    # cat('No transform done. \n')
  }

  return(predictions)

}


#' @title Predict TF occupancy using posterior mean of regression coefficients
#' @description Predict TF occupancy using posterior mean of regression
#' coefficients trained from the TOP model
#'
#' @param data A data frame. Columns are motif score and DNase (or ATAC) bins.
#' Rows are candidate sites.
#' @param mean_coef A numeric vector. The posterior mean of trained regression
#' coefficients, including the intercept and coefficients for motif score and
#' DNase features.
#' length(mean_coef) should be equal to 1+ncol(data).
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Options: asinh, log2, log, none.
#'
#' @return The function returns a vector of predicted TF occupancy.
#'
#' @export
predict_TOP_mean_coef <- function(data,
                                  mean_coef,
                                  transform = c('asinh', 'log2', 'log', 'none')){

  transform <- match.arg(transform)

  if((ncol(data)+1) != length(mean_coef)){
    stop('The number of coefficients should be equal to
         the number of data columns + 1 (intercept)!')
  }

  coefficients <- as.matrix(mean_coef, ncol = 1)

  data <- as.matrix(data.frame(intercept = 1, data, check.names = F))

  predictions <- as.numeric(data %*% coefficients)

  if(transform == 'asinh'){
    predictions <- sinh(predictions)
  }else if (transform == 'log2'){ # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'log'){ # log(y+1)
    predictions <- exp(predictions) - 1
  }else if (transform == 'none'){ # no transform
    # cat('Do not transform predictions. \n')
  }

  return(predictions)

}

#' @title Predict TF binding by TOP logistic model
#' with posterior mean of regression coefficients
#' @description Predict TF binding using posterior mean of the regression
#' coefficients trained from TOP logistic model
#'
#' @param data A data frame. Columns are motif score and DNase (or ATAC) bins.
#' Rows are candidate sites.
#' @param mean_coef A numeric vector. The posterior mean of trained regression
#' coefficients, including the intercept and coefficients for motif score and
#' DNase features.
#' length(mean_coef) should be equal to 1+ncol(data).
#'
#' @return The function returns a vector of predicted TF binding probabilities
#'
#' @export
predict_TOP_logistic_mean_coef <- function(data, mean_coef){

  if((ncol(data)+1) != length(mean_coef)){
    stop('The number of coefficients should be equal to
         the number of data columns + 1 (intercept)!')
  }

  coefficients <- as.matrix(mean_coef, ncol = 1)

  data <- as.matrix(data.frame(intercept = 1, data, check.names = F))

  mu <- as.numeric(data %*% coefficients)

  p <- 1 / (1 + exp(-mu)) # inverse logit

  return(p)

}

