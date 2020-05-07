
#' @title Predict TF occupancy using posterior samples of regression coefficients
#' @description Predict TF occupancy using posterior samples of regression
#' coefficients trained from the TOP model
#'
#' @param data A data frame or matrix. Columns are motif score and DNase features.
#' Rows are candidate sites.
#' @param alpha_samples posterior samples of alpha
#' @param beta_samples posterior samples of beta
#' @param tau_samples posterior samples of tau
#' @param sample logicals. If TRUE, samples from posterior predictions and then
#' take the mean of posterior prediction samples.
#' @param average_parameters logicals. If TRUE, uses the posterior mean of
#' regression coefficients to make predictions.
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Default: transform = "asinh".
#'
#' @return The function returns a vector of predicted TF occupancy.
#' @export
#'
#' @examples
#' predicted <- predict_norm(data, alpha_samples, beta_samples, tau_samples,
#' sample = TRUE, average_parameters = FALSE, transform = 'asinh')
#'
predict_norm <- function(data, alpha_samples, beta_samples, tau_samples,
                         sample = TRUE, average_parameters = FALSE, transform = 'asinh'){

  # the mean of the posterior norm is x %*% beta + alpha

  coefficients <- t(cbind(alpha_samples, beta_samples))

  data <- as.matrix(data.frame(intercept = 1, data, check.names = F))

  if(average_parameters){
    coefficients <- apply(coefficients, 1, mean)
    means <- data %*% coefficients
    predictions <- as.numeric(means)
  }else{
    means <- data %*% coefficients

    if(sample){
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
  }else if (transform == 'log2'){ # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'log'){ # log(y+1)
    predictions <- exp(predictions) - 1
  }else if (transform == ''){ # no transform
    cat('Do not transform predictions. \n')
  }else{
    warning('No transform done. Please check the transform method! \n')
  }

  return(predictions)
}


#' @title Predict TF occupancy using posterior mean of regression coefficients
#' @description Predict TF occupancy using posterior mean of regression
#' coefficients trained from the TOP model
#'
#' @param data A data frame or matrix. Columns are motif score and DNase features.
#' Rows are candidate sites.
#' @param coef_mean A numeric vector. The posterior mean of trained regression
#' coefficients, including the intercept and coefficients for motif score and
#' DNase features.
#' length(coef_mean) should be equal to 1+ncol(data).
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Default: transform = "asinh".
#'
#' @return The function returns a vector of predicted TF occupancy.
#' @export
#'
#' @examples
#' predicted <- predict_coef_BH_mean(data, coef_mean, transform = 'asinh')
#'
predict_coef_BH_mean <- function(data, coef_mean, transform = 'asinh'){

  if((ncol(data)+1) != length(coef_mean)){
    stop('The number of coefficients is not equal to the number of data columns + 1!')
  }

  # the mean of the posterior norm is alpha + X %*% beta
  coefficients <- as.matrix(coef_mean, ncol = 1)

  data <- as.matrix(data.frame(intercept = 1, data, check.names = F))

  predictions <- as.numeric(data %*% coefficients)

  if(transform == 'asinh'){
    predictions <- sinh(predictions)
  }else if (transform == 'log2'){ # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'log'){ # log(y+1)
    predictions <- exp(predictions) - 1
  }else if (transform == ''){ # no transform
    cat('Do not transform predictions. \n')
  }else{
    warning('No transform done. Please check the transform method! \n')
  }

  return(predictions)

}

