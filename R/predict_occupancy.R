
#' @title Predict quantitative TF occupancy or TF binding probability
#' @description Predict quantitative TF occupancy or TF binding probability
#' using TOP model trained from ChIP-seq read counts or binary labels.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param tf_name specifies TF name.
#' @param cell_type specifies the cell type.
#' @param TOP_coef A list containing the posterior mean of TOP regression coefficients.
#' @param level Specific the TOP model hierarchy level to use.
#' Options: 'best', 'bottom', 'middle', or 'top'.
#' Default: 'best' -- choose the best model (lowest level of the hierarchy available):
#' If the TF motif and cell type is available in the training data,
#' then use the bottom level (TF- and cell-type-specific model).
#' otherwise, if TF motif (but not cell type) is available in the training data,
#' choose the middle level (TF-specific model) of that TF motif;
#' otherwise, use the top level TF-generic model.
#' @param logistic.model Logical; if TRUE, use the logistic version of TOP model.
#' @param transform Type of transformation for ChIP counts.
#' Possible values are "asinh", "log2", "sqrt", and "none" (no transformation).
#' Only needed when logistic.model is FALSE.
#'
#' @export
#'
predict_TOP <- function(data,
                        tf_name,
                        cell_type,
                        TOP_coef,
                        level = c('best', 'bottom', 'middle', 'top'),
                        logistic.model = FALSE,
                        transform = c('asinh', 'log2', 'log', 'none')){

  level <- match.arg(level)

  if(missing(tf_name)){
    tf_name = NA
  }
  if(missing(cell_type)){
    cell_type = NA
  }

  selected.model <- select_model_coef_level(tf_name, cell_type, TOP_coef, level)

  if(logistic.model == FALSE) {
    transform <- match.arg(transform)
    predicted <- predict_TOP_mean_coef(data, selected.model$coef, transform = transform)
  }else if (logistic.model == TRUE){
    predicted <- predict_TOP_logistic_mean_coef(data, selected.model$coef)
  }

  res <- list(level = selected.model$level,
              model = selected.model$model,
              coef = selected.model$coef,
              data = data,
              predicted = predicted)
  return(res)

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

  cat('Predicting TF occupancy using TOP occupancy model...\n')
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

  cat('Predicting TF occupancy using TOP occupancy model with posterior samples...\n')
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

  cat('Predicting TF binding probability using TOP logistic model...\n')

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

#' @title Select PWM and DNase (or ATAC) bin features from the input data frame
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


#' @title Select regression coefficients by the TOP hierarchy level
#'
#' @param tf_name specifies TF name.
#' @param cell_type specifies the cell type.
#' @param TOP_mean_coef Trained TOP posterior mean regression coefficients.
#' @param level Specific the TOP model hierarchy level to use.
#' Options: 'best', 'bottom', 'middle', or 'top'.
#' Default: 'best' -- choose the best model (lowest level of the hierarchy available):
#' If the TF motif and cell type is available in the training data,
#' then use the bottom level (TF- and cell-type-specific model).
#' otherwise, if TF motif (but not cell type) is available in the training data,
#' choose the middle level (TF-specific model) of that TF motif;
#' otherwise, use the top level TF-generic model.
#' @export
#'
select_model_coef_level <- function(tf_name,
                                    cell_type,
                                    TOP_mean_coef,
                                    level = c('best', 'bottom', 'middle', 'top')) {

  level <- match.arg(level)

  bottom_level_mean_coef <- TOP_mean_coef$bottom
  middle_level_mean_coef <- TOP_mean_coef$middle
  top_level_mean_coef <- TOP_mean_coef$top

  if(level == 'best'){

    ## load model, using lower level model if available
    tf_cell_name <- paste(tf_name, cell_type, sep = '.')
    if (tf_cell_name %in% rownames(bottom_level_mean_coef)) {
      model_coef <- bottom_level_mean_coef[tf_cell_name, ]
      model_level <- 'bottom'
      model <- tf_cell_name
    } else if (tf_name %in% rownames(middle_level_mean_coef)) {
      model_coef <- middle_level_mean_coef[tf_name, ]
      model_level <- 'middle'
      model <- tf_name
    } else {
      model_coef <- top_level_mean_coef
      model_level <- 'top'
      model <- 'TF-generic'
    }
    cat(model, model_level, 'level model selected. \n')

  }else if (level == 'bottom'){

    tf_cell_name <- paste(tf_name, cell_type, sep = '.')

    if (tf_cell_name %in% rownames(bottom_level_mean_coef)) {
      model_coef <- bottom_level_mean_coef[tf_cell_name, ]
      model_level <- 'bottom'
      model <- tf_cell_name
      cat(model, model_level, 'level model selected. \n')
    } else{
      model_coef <- NA
      model_level <- 'bottom'
      model <- NA
      cat(model_level, 'level model is not available! \n')
    }

  }else if (level == 'middle'){

    if (tf_name %in% rownames(middle_level_mean_coef)) {
      model_coef <- middle_level_mean_coef[tf_name, ]
      model_level <- 'middle'
      model <- tf_name
      cat(model, model_level, 'level model selected. \n')
    } else{
      model_coef <- NA
      model_level <- 'middle'
      cat(model_level, 'level model is not available! \n')
    }

  }else if(level == 'top'){

    cat('Choose top level model. \n')
    model_coef <- top_level_mean_coef
    model_level <- 'top'
    model <- 'TF-generic'
    cat(model, model_level, 'level model selected. \n')
  }

  return(list(level = model_level,
              model = model,
              coef = model_coef))

}
