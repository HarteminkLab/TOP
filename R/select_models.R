#' @title Select TOP model hierarchy level
#' @description If the TF motif is available in the training data and
#' model_level='best', then use the TF-specific model of that TF motif;
#' otherwise, use the top level TF-generic model.
#'
#' @param pwm_id Motif PWM ID
#' @param cell_type Cell type
#' @param tf_pwm_training_list List of TF motifs with available coefficients trained.
#' @param model_coefficients Trained TOP model regression coefficients
#' @param model_level Specific the TOP model hierarchy level:
#' 'top', 'bottom', 'middle', or 'best'.
#' Default: 'best', the lowest hierarchy level of the model will be used.
#'
select_model_coef_mean <- function(pwm_id,
                                   cell_type,
                                   tf_pwm_training_list,
                                   model_coefficients,
                                   model_level = c('best','bottom', 'middle', 'top') ) {

  model_level <- match.arg(model_level)

  coef_TFs_bottom_mean.m <- model_coefficients$bottom
  coef_TFs_middle_mean.m <- model_coefficients$middle
  coef_TFs_top_mean.m <- model_coefficients$top

  pwm_stable_id <- gsub('[.].+', '', pwm_id)

  if (pwm_stable_id %in% tf_pwm_training_list$pwm_stable_id && model_level == 'best') {
    cat('TF model available. \n')
    tf_model <- tf_pwm_training_list[tf_pwm_training_list$pwm_stable_id == pwm_stable_id, 'tf_name']
  } else {
    cat('No TF-specific model available. Use TF generic model. ')
    tf_model <- 'TF-generic'
  }

  cat('Use', tf_model, 'model in', cell_type, '\n')

  ## load model, using lower level model if available
  tf_cell_name <- paste(tf_model, cell_type, sep = '.')

  if (tf_cell_name %in% rownames(coef_TFs_bottom_mean.m)) {
    cat('Use bottom level model. \n')
    coef_mean_tf <- coef_TFs_bottom_mean.m[tf_cell_name, ]
    model_type <- 'Bottom'

  } else if (tf_model %in% rownames(coef_TFs_middle_mean.m)) {
    cat('Use middle level model. \n')
    coef_mean_tf <- coef_TFs_middle_mean.m[tf_model, ]
    model_type <- 'Middle'

  } else {
    cat('Use top level model. \n')
    coef_mean_tf <- coef_TFs_top_mean.m
    model_type <- 'Top'

  }

  return(list(tf_model = tf_model, model_type = model_type, coef_mean_tf = coef_mean_tf))

}
