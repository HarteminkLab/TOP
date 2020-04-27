
#' @title Predict TF occupancy using posterior mean of regression coefficients
#' @description Predict TF occupancy using posterior mean of regression
#' coefficients trained using the TOP model
#'
#' @param data A data frame or matrix. Columns are motif score and DNase features.
#' Rows are candidate sites.
#' @param coef_mean A numeric vector. The posterior mean of trained regression
#' coefficients, including the intercept and coefficients for motif score and
#' DNase features.
#' length(coef_mean) should be equal to 1+ncol(data).
#' @param transform Method used to transform ChIP-seq occupancy when training
#' the TOP model.
#' Default: transform = "asinh".
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
