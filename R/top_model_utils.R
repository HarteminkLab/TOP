#' @title Select TOP posterior samples
#' @description Select posterior samples and perform thinning and sampling
#' of the posterior samples if needed.
#' @param TOP_samples TOP posterior samples.
#' @param thin Thinning rate of extract the posterior samples,
#' must be a positive integer (default = 1, no thinning performed).
#' @param n_samples Keep n_samples posterior samples (randomly choose),
#' when the number of posterior samples is greater than \code{n_samples}.
#' @return A data frame of posterior samples.
#' @export
#'
select_TOP_samples <- function(TOP_samples, thin = 1, n_samples = 1000) {

  # If samples are from multiple MCMC chains, combine samples from the chains
  if ( is.list(TOP_samples) ){
    TOP_samples <- as.data.frame(do.call(rbind, TOP_samples))
  }

  if (thin > 1) {
    TOP_samples <- TOP_samples[seq(from = 1, to = nrow(TOP_samples), by = thin), ]
  }

  if (nrow(TOP_samples) > n_samples) {
    TOP_samples <- TOP_samples[seq(from = 1, to = nrow(TOP_samples), length.out = n_samples),]
  }

  return(TOP_samples)

}

#' @title Combine and take the average of TOP posterior samples from all partitions
#' @description Combine and take the average of TOP posterior samples
#' from all partitions. Use \code{select_TOP_samples} to select
#' posterior samples from each partition.
#' @param TOP_samples_files Files of TOP posterior samples from all partitions.
#' @param thin thinning rate of extract the posterior samples,
#' must be a positive integer (default = 1, no thinning performed).
#' @param n_samples Keep n_samples posterior samples (randomly choose),
#' when the number of posterior samples is greater than \code{n_samples}.
#' @return A data frame of combined and averaged posterior samples.
#' @export
#' @examples
#' TOP_samples <- combine_TOP_samples(TOP_samples_files, n_samples = 1000)
combine_TOP_samples <- function(TOP_samples_files,
                                thin = 1,
                                n_samples = 1000) {
  if(any(!file.exists(TOP_samples_files))){
    stop('Files not available:\n', TOP_samples_files[!file.exists(TOP_samples_files)])
  }
  N <- length(all_TOP_samples)
  cat('Loading samples from partition 1 ...\n')
  TOP_samples <- readRDS(TOP_samples_files[1])
  combined_TOP_samples <- select_TOP_samples(all_TOP_samples[[1]], thin, n_samples)
  if(N > 1) {
    for(i in 2:N) {
      cat('Loading samples from partition', i, '...\n')
      TOP_samples <- readRDS(TOP_samples_files[i])
      combined_TOP_samples <- combined_TOP_samples + select_TOP_samples(TOP_samples, thin, n_samples)
    }
    combined_TOP_samples <- combined_TOP_samples / N
  }
  return(combined_TOP_samples)
}

#' @title Extract alpha and beta coefficients from TOP posterior samples
#'
#' @param TOP_samples TOP samples combined from all partitions using
#' \code{combine_TOP_samples}.
#' @param tf_cell_combos A table with the indices and names of TF and cell type
#' combinations from the assembled training data. If missing,
#' will extract from assembled_training_data.
#' @param assembled_training_data Assembled training data as in the
#' \code{assemble_training_data} function.
#' @param tf_name TF name.
#' @param cell_type Cell type.
#' @param n_bins Number of DNase or ATAC bins in TOP model (default = 5).
#' @param level The level in the TOP model (options: bottom, middle, or top),
#' 'bottom' level: TF- and cell-type- specific,
#' 'middle' level: TF-specific, cell-type generic,
#' 'top' level: TF-generic
#' @return A data frame of posterior samples for TOP's alpha and beta coefficients.
#' @export
#'
extract_TOP_coef_samples <- function(TOP_samples,
                                     tf_cell_combos,
                                     assembled_training_data,
                                     tf_name,
                                     cell_type,
                                     n_bins = 5,
                                     level = c('bottom', 'middle', 'top')){

  level <- match.arg(level)

  if(missing(tf_cell_combos)){
    if(!missing(assembled_training_data)){
      tf_cell_combos <- extract_tf_cell_combos(assembled_training_data)
    }else{
      stop('Please provide tf_cell_combos or assembled_training_data!')
    }
  }

  TOP_samples <- as.data.frame(TOP_samples)

  if(level == 'bottom' && !missing(tf_name) && !missing(cell_type)){
    tf_cell_combo <- tf_cell_combos[tf_cell_combos$tf_name == tf_name & tf_cell_combos$cell_type == cell_type,]
    tf_id <- tf_cell_combo$tf_id[1]
    cell_id <- tf_cell_combo$cell_id[1]

    # alpha's and beta's are TF- and cell-type-specific model coefficients
    # specific for each TF and cell type combination
    data_id <- sprintf('[%s,%s]', tf_id, cell_id)
    coef_names <- c(paste0('alpha', data_id), paste0('beta', 1:(1+n_bins), data_id))
    TOP_coef_samples <- TOP_samples[, coef_names]

  }else if( level == 'middle' && !missing(tf_name) ){
    tf_id <- tf_cell_combos$tf_id[tf_cell_combos$tf_name == tf_name][1]

    # Alpha and Beta's are TF specific, cell-type-generic model coefficients
    data_id <- sprintf('[%s]', tf_id)
    coef_names <- c(paste0('Alpha', data_id), paste0('Beta', 1:(1+n_bins), data_id))
    TOP_coef_samples <- TOP_samples[, coef_names]

  }else if(level == 'top'){
    # A and B's are TF-generic model coefficients
    coef_names <- c('A', paste0('B', 1:(1+n_bins)))
    TOP_coef_samples <- TOP_samples[, coef_names]

  }

  return(TOP_coef_samples)

}

#' @title Extract the posterior mean coefficients for each level of TOP model
#' @description Extract alpha and beta coefficients from TOP posterior samples,
#' and return the mean of coefficients at each level of TOP model.
#' @param TOP_samples TOP samples combined from all partitions using
#' \code{combine_TOP_samples}.
#' @param tf_cell_combos A table with the indices and names of TF and cell type
#' combinations from the assembled training data. If missing,
#' will extract from assembled_training_data.
#' @param assembled_training_data Assembled training data as in the
#' \code{assemble_training_data} function.
#' @param n_bins Number of DNase or ATAC bins in TOP model (default = 5)
#' @return A list of posterior mean coefficients at each level of TOP model.
#' @export
#'
extract_TOP_mean_coef <- function(TOP_samples,
                                  tf_cell_combos,
                                  assembled_training_data,
                                  n_bins = 5){

  # Extract the indices and names of TF and cell type combos
  if(missing(tf_cell_combos)){
    if(!missing(assembled_training_data)){
      tf_cell_combos <- extract_tf_cell_combos(assembled_training_data)
    }else{
      stop('Please provide tf_cell_combos or assembled_training_data!')
    }
  }

  # Extract bottom level posterior mean coefficients for all TF x cell type combos
  bottom_level_mean_coef <- matrix(NA, nrow = nrow(tf_cell_combos), ncol = 2+n_bins)
  for( i in 1:nrow(tf_cell_combos)){
    bottom_level_coef_samples <- extract_TOP_coef_samples(TOP_samples,
                                                          tf_cell_combos,
                                                          tf_name = tf_cell_combos[i, 'tf_name'],
                                                          cell_type = tf_cell_combos[i, 'cell_type'],
                                                          n_bins = n_bins,
                                                          level = 'bottom')
    bottom_level_mean_coef[i, ] <- colMeans(bottom_level_coef_samples)
  }
  rownames(bottom_level_mean_coef) <- paste(tf_cell_combos$tf_name, tf_cell_combos$cell_type, sep = '.')
  colnames(bottom_level_mean_coef) <- c('Intercept', 'PWM', paste0('Bin', 1:n_bins))

  # Extract middle level posterior mean coefficients for all TFs
  tf_name_list <- unique(as.character(tf_cell_combos$tf_name))
  middle_level_mean_coef <- matrix(NA, nrow = length(tf_name_list), ncol = 2+n_bins)
  for( i in 1:length(tf_name_list)){
    middle_level_coef_samples <- extract_TOP_coef_samples(TOP_samples,
                                                          tf_cell_combos,
                                                          tf_name = tf_name_list[i],
                                                          n_bins = n_bins,
                                                          level = 'middle')
    middle_level_mean_coef[i, ] <- colMeans(middle_level_coef_samples)
  }
  rownames(middle_level_mean_coef) <- tf_name_list
  colnames(middle_level_mean_coef) <- c('Intercept', 'PWM', paste0('Bin', 1:n_bins))

  # Extract top level posterior mean coefficients
  top_level_coef_samples <- extract_TOP_coef_samples(TOP_samples,
                                                     tf_cell_combos,
                                                     n_bins = n_bins,
                                                     level = 'top')
  top_level_mean_coef <- colMeans(top_level_coef_samples)
  names(top_level_mean_coef) <- c('Intercept', 'PWM', paste0('Bin', 1:n_bins))

  # Combine posterior mean coefficients for all three levels
  TOP_mean_coef <- list(top = top_level_mean_coef,
                        middle = middle_level_mean_coef,
                        bottom = bottom_level_mean_coef)

  return(TOP_mean_coef)
}

#' @title Extract a table listing the indices and names of TF and cell type
#' combinations from the assembled training data
#' @description Extract a table listing the indices and names of TF x cell type
#' combinations from the assembled training data.
#' This will be used to extract regression coefficients from the TOP fit result.
#' @param assembled_training_data Assembled training data as in the
#' \code{assemble_training_data} function
#'
#' @return a data frame of the indices and names of TF x cell type combinations.
#' @export
extract_tf_cell_combos <- function(assembled_training_data){
  cat('Extract the table of TF x cell combinations from assembled_training_data...\n')
  if(inherits(assembled_training_data, "list"))
    assembled_training_data <- assembled_training_data[[1]]

  tf_cell_combos <- unique(assembled_training_data[, c('tf_id', 'cell_id', 'tf_name', 'cell_type')])
  return(tf_cell_combos)
}
