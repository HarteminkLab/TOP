
#' @title Predict quantitative TF occupancy or TF binding probability
#' @description Predict quantitative TF occupancy or TF binding probability
#' using TOP model trained from ChIP-seq read counts or binary labels.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param TOP_coef A list containing the posterior mean of TOP regression coefficients.
#' @param tf_name specifies the TF name to make predictions.
#' It will find the model parameters trained for this TF.
#' This is not needed (not used) when \code{level = 'top'}.
#' @param cell_type specifies the cell type to make predictions.
#' It will find the model parameters trained for this cell type.
#' This is not needed (not used) when \code{level = 'middle'} or \code{level = 'top'}.
#' @param level Specific the TOP model hierarchy level to use.
#' Options: \sQuote{best}, \sQuote{bottom}, \sQuote{middle}, or \sQuote{top}.
#' When \code{level = 'best'}, use the best (lowest available) level of the
#' hierarchy for the TF x cell type combination.
#' When \code{level = 'bottom'}, use the bottom level (TF- and cell-type-specific model),
#' if the TF motif and cell type is available in the training data.
#' When \code{level = 'middle'}, use the middle level (TF-specific model) of that TF.
#' When \code{level = 'top'}, use the top level TF-generic model.
#' @param logistic_model Logical. Whether to use the logistic version of TOP model.
#' If \code{logistic_model = TRUE},
#' use the logistic version of TOP model to predict TF binding probability.
#' If \code{logistic_model = FALSE}, use the quantitative occupancy model (default).
#' @param transform Type of transformation performed for ChIP-seq read counts
#' when preparing the input training data.
#' Options are: \sQuote{asinh}(asinh transformation),
#' \sQuote{log2} (log2 transformation),
#' \sQuote{sqrt} (square root transformation),
#' and \sQuote{none}(no transformation).
#' This only applies when \code{logistic_model = FALSE}.
#' @return Returns a list with the following elements,
#'
#' \item{model}{TOP model name.}
#' \item{level}{selected hierarchy level.}
#' \item{coef}{posterior mean of regression coefficients.}
#' \item{predictions}{a data frame with the data and predicted values.}
#' @examples
#' # 'data' is a data frame of input data,
#' # with columns of PWM scores and five DNase (or ATAC) bins.
#' # 'TOP_coef' is pretrained posterior mean of TOP regression coefficients.
#'
#' # Predict CTCF occupancy in K562 using the quantitative occupancy model:
#'
#' # Predict using the 'bottom' level model
#' result <- predict_TOP(data, TOP_coef,
#'                       tf_name = 'CTCF', cell_type = 'K562',
#'                       level = 'bottom',
#'                       logistic_model = FALSE,
#'                       transform = 'asinh') # We used 'asinh' transformation on the ChIP data when training the model.
#'
#' # Predict using the 'best' model
#' # Since CTCF in K562 cell type is included in training,
#' # the 'best' model is the 'bottom' level model.
#' result <- predict_TOP(data, TOP_coef,
#'                       tf_name = 'CTCF', cell_type = 'K562', level = 'best',
#'                       logistic_model = FALSE, transform = 'asinh')
#'
#' # We can use the 'middle' model to predict CTCF in K562
#' # or other cell types or conditions
#' result <- predict_TOP(data, TOP_coef,
#'                       tf_name = 'CTCF', level = 'middle',
#'                       logistic_model = FALSE, transform = 'asinh')
#'
#' # Predict CTCF binding probability using the logistic version of the model:
#' # No need to set the argument for 'transform' for the logistic model.
#'
#' # Predict using the 'bottom' level model
#' result <- predict_TOP(data, TOP_coef,
#'                       tf_name = 'CTCF', cell_type = 'K562',
#'                       level = 'best',
#'                       logistic_model = TRUE)
#'
#' # Predict using the 'middle' level model
#' result <- predict_TOP(data, TOP_coef,
#'                       tf_name = 'CTCF', level = 'middle',
#'                       logistic_model = TRUE)
#' @export
#'
predict_TOP <- function(data,
                        TOP_coef,
                        tf_name,
                        cell_type,
                        level = c('best', 'bottom', 'middle', 'top'),
                        logistic_model = FALSE,
                        transform = c('asinh', 'log2', 'log', 'none')){

  level <- match.arg(level)

  if (missing(TOP_coef)){
    stop("Please provide TOP regression coefficients in 'TOP_coef'.")
  }

  selected_model <- select_model_coef_level(TOP_coef, tf_name, cell_type, level)
  if (anyNA(selected_model$coef)){
    stop(sprintf("%s level cefficients are not available.\n", level))
  }

  if (logistic_model == FALSE ) {
    transform <- match.arg(transform)
    predictions <- predict_TOP_mean_coef(data, selected_model$coef, transform = transform)
  }else if (logistic_model == TRUE){
    predictions <- predict_TOP_logistic_mean_coef(data, selected_model$coef)
  }

  return(list(model = selected_model$model,
              level = selected_model$level,
              coef = selected_model$coef,
              predictions = predictions))

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
#' @return A data frame of input data and predicted TF occupancy (posterior mean).
#'
#' @export
#' @examples
#' predictions <- predict_TOP_mean_coef(data, # data frame of input data.
#'                                      mean_coef, # (selected) pretrained posterior mean of TOP regression coefficients
#'                                      transform = 'asinh') # We used 'asinh' transformation on ChIP counts when training the model
#'
predict_TOP_mean_coef <- function(data,
                                  mean_coef,
                                  transform = c('asinh', 'log2', 'log', 'none')){

  transform <- match.arg(transform)

  features <- select_features(data)

  if((ncol(features)+1) != length(mean_coef)){
    stop('The number of coefficients not equal to
         the number of features + intercept! Check input data!')
  }

  cat('Predicting TF occupancy using TOP occupancy model...\n')

  coefficients <- as.matrix(mean_coef, ncol = 1)

  data_matrix <- as.matrix(data.frame(intercept = 1, features, check.names = FALSE))

  predictions <- as.numeric(data_matrix %*% coefficients)

  # transform back to the original scale
  if(transform == 'asinh'){
    predictions <- sinh(predictions)
  }else if (transform == 'log2'){
    # log2(y+1)
    predictions <- 2^predictions - 1
  }else if (transform == 'sqrt'){
    # sqrt(y)
    predictions <- predictions^2
  }

  return(data.frame(data, predicted = predictions))

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
#' @return A data frame of input data and predicted TF binding probability.
#'
#' @export
#' @examples
#' predictions <- predict_TOP_logistic_mean_coef(data, # data frame of input data.
#'                                               mean_coef)  # (selected) pretrained posterior mean of TOP regression coefficients
#'
predict_TOP_logistic_mean_coef <- function(data, mean_coef){

  features <- select_features(data)

  if((ncol(features)+1) != length(mean_coef)){
    stop('The number of coefficients not equal to
         the number of features + intercept! Check input data!')
  }

  cat('Predicting TF binding probability using TOP logistic model...\n')

  coefficients <- as.matrix(mean_coef, ncol = 1)

  data_matrix <- as.matrix(data.frame(intercept = 1, features, check.names = FALSE))

  mu <- as.numeric(data_matrix %*% coefficients)

  p <- 1 / (1 + exp(-mu)) # inverse logit

  return(data.frame(data, predicted = p))


}


#' @title Predict TF occupancy using posterior samples of regression coefficients
#' @description Predict TF occupancy using posterior samples of TOP regression
#' coefficients.
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param coef_samples TOP posterior samples.
#' @param use_posterior_mean Logical; if TRUE, uses the posterior mean of
#' regression coefficients to make predictions.
#' @param sample_predictions Logical; if TRUE, sample from posterior predictions
#' and then take the mean of posterior prediction samples.
#' @param transform Method used to transform ChIP-seq counts when training
#' the TOP model. Options: asinh, log2, sqrt, none.
#'
#' @return A data frame of input data and predicted TF occupancy (posterior mean).
#' @export
#' @examples
#' predictions <- predict_TOP_samples(data,
#'                                    coef_samples = TOP_samples,
#'                                    use_posterior_mean = FALSE,
#'                                    sample_predictions = TRUE,
#'                                    transform = "asinh")
#'
predict_TOP_samples <- function(data,
                                coef_samples,
                                use_posterior_mean = FALSE,
                                sample_predictions = TRUE,
                                transform = c('asinh', 'log2', 'sqrt', 'none')){

  transform <- match.arg(transform)

  features <- select_features(data)

  cat('Predicting TF occupancy using TOP occupancy model...\n')
  data_matrix <- as.matrix(data.frame(intercept = 1, features, check.names = FALSE))

  alpha_samples <- coef_samples$alpha_samples
  beta_samples  <- coef_samples$beta_samples
  tau_samples   <- coef_samples$tau_samples
  coefficients  <- t(cbind(alpha_samples, beta_samples))

  if(use_posterior_mean){
    coefficients <- apply(coefficients, 1, mean)
    means <- data_matrix %*% coefficients
    predictions <- as.numeric(means)
  }else{
    means <- data_matrix %*% coefficients

    if(sample_predictions){
      sds <- 1 / sqrt(tau_samples)
      n_data <- nrow(data_matrix)
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
  }

  return(data.frame(data, predicted = predictions))

}


#' @title Select regression coefficients by the TOP hierarchy level
#'
#' @param TOP_mean_coef Trained TOP posterior mean regression coefficients.
#' @param tf_name specifies TF name.
#' @param cell_type specifies the cell type.
#' @param level Specific the TOP model hierarchy level to use.
#' Options: 'best', 'bottom', 'middle', or 'top'.
#' Default: 'best' -- choose the best model (lowest level of the hierarchy available):
#' If the TF motif and cell type is available in the training data,
#' then use the bottom level (TF- and cell-type-specific model).
#' otherwise, if TF motif (but not cell type) is available in the training data,
#' choose the middle level (TF-specific model) of that TF motif;
#' otherwise, use the top level TF-generic model.
#' @return A list containing the results, including selected hierarchy level,
#' TOP model, regression coefficients (posterior mean).
#' @export
select_model_coef_level <- function(TOP_mean_coef,
                                    tf_name,
                                    cell_type,
                                    level = c('best', 'bottom', 'middle', 'top')) {

  level <- match.arg(level)

  if( missing(tf_name) && level %in% c('best', 'bottom', 'middle') ){
    stop(sprintf("Please specify 'tf_name' when 'level = %s'.", level))
  }

  if( missing(cell_type) && level %in% c('best', 'bottom')){
    stop(sprintf("Please specify 'cell_type' when 'level = %s'.", level))
  }

  # convert TF name to upper case (as we use upper case for TF names in training)
  if ( level %in% c('best', 'bottom', 'middle') ){
    tf_name <- base::toupper(tf_name)
  }

  bottom_level_mean_coef <- TOP_mean_coef$bottom
  middle_level_mean_coef <- TOP_mean_coef$middle
  top_level_mean_coef <- TOP_mean_coef$top

  if(level == 'best'){
    ## load model, using lower level model if available
    tf_cell_name <- paste(tf_name, cell_type, sep = '.')
    if (tf_cell_name %in% rownames(bottom_level_mean_coef)) {
      model_coef <- bottom_level_mean_coef[tf_cell_name, ]
      model_level <- 'bottom'
      model_name <- tf_cell_name
      cat(sprintf('Use the bottom level model for %s in %s cell type.\n', tf_name, cell_type))
    } else if (tf_name %in% rownames(middle_level_mean_coef)) {
      model_coef <- middle_level_mean_coef[tf_name, ]
      model_level <- 'middle'
      model_name <- tf_name
      cat(sprintf('Use the middle level model for %s.\n', tf_name))
    } else {
      model_coef <- top_level_mean_coef
      model_level <- 'top'
      model_name <- 'TF-generic'
      cat('Use the top level model.\n')
    }
  }else if (level == 'bottom'){
    tf_cell_name <- paste(tf_name, cell_type, sep = '.')
    if (tf_cell_name %in% rownames(bottom_level_mean_coef)) {
      model_coef <- bottom_level_mean_coef[tf_cell_name, ]
      model_level <- 'bottom'
      model_name <- tf_cell_name
      cat(sprintf('Use the bottom level model for %s in %s cell type.\n', tf_name, cell_type))
    } else{
      model_coef <- NA
      model_level <- 'bottom'
      model_name <- NA
      cat(model_level, 'level model is not available! \n')
    }
  }else if (level == 'middle'){
    if (tf_name %in% rownames(middle_level_mean_coef)) {
      model_coef <- middle_level_mean_coef[tf_name, ]
      model_level <- 'middle'
      model_name <- tf_name
      cat(sprintf('Use the middle level model for %s.\n', tf_name))
    } else{
      model_coef <- NA
      model_level <- 'middle'
      cat(model_level, 'level model is not available! \n')
    }
  }else if(level == 'top'){
    model_coef <- top_level_mean_coef
    model_level <- 'top'
    model_name <- 'TF-generic'
    cat('Use the top level model.\n')
  }

  return(list(level = model_level,
              model = model_name,
              coef = model_coef))

}


#' @title Select PWM and DNase (or ATAC) bin features from the input data
#'
#' @param data A data frame containing motif PWM score and DNase (or ATAC) bins.
#' @param pwm_col name (prefix) of the PMW score column.
#' @param bin_col name (prefix) of the DNase or ATAC bin columns.
#' @param quiet Logical. If TRUE, disable messages.
select_features <- function(data, pwm_col = 'pwm', bin_col = 'bin', quiet=FALSE){
  data <- as.data.frame(data)
  pwm_col <- grep(pwm_col, colnames(data), ignore.case = TRUE, value = TRUE)
  bin_cols <- grep(bin_col, colnames(data), ignore.case = TRUE, value = TRUE)
  features_cols <- c(pwm_col, bin_cols)
  if(!quiet){
    cat('Select features:', features_cols, '\n')
  }
  features <- data[, features_cols]
  return(features)
}

