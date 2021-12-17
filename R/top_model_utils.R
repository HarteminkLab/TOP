
#' Load TOP posterior samples
#'
#' @param TOP_samples_file file names of the posterior samples from all partitions
#' @param thin thinning rate of extract the posterior samples,
#' must be a positive integer.
#' @param n.samples Randomly choose n.samples posterior samples,
#' if the number of posterior samples is greater than n.samples
#'
#' @return
#' @export
#'
load_TOP_samples <- function(TOP_samples_file, thin = 1, n.samples = 1000) {

  if (!file.exists(TOP_samples_file)) {
    stop(paste("TOP sample file: ", TOP_samples_file, "is not available!"))
  }
  TOP_samples <- readRDS(TOP_samples_file)

  # merge multiple MCMC chains
  if ( is.list(TOP_samples) ){
    TOP_samples <- as.data.frame(do.call(rbind, TOP_samples))
  }

  if (thin > 1) {
    TOP_samples <- TOP_samples[seq(from = 1, to = nrow(TOP_samples), by = thin), ]
  } else if (nrow(TOP_samples) > n.samples) {
    TOP_samples <- TOP_samples[seq(from = 1, to = nrow(TOP_samples), length.out = n.samples),]
  }

  return(TOP_samples)

}

#' Combine and take the average of TOP posterior samples from all partitions
#'
#' @param TOP_samples_files file names of the posterior samples from all partitions
#' @param thin thinning rate of extract the posterior samples,
#' must be a positive integer.
#' @param n.samples Randomly choose n.samples posterior samples,
#' if the number of posterior samples is greater than n.samples
#'
#' @return
#' @export
#'
combine_TOP_samples <- function(TOP_samples_files, thin = 1, n.samples = 1000) {

  N <- length(TOP_samples_files)
  cat('Combining', N, 'sample partitions ... \n')
  TOP_samples <- load_TOP_samples(TOP_samples_files[1], thin, n.samples)

  if(N > 1) {
    for(i in 2:N) {
      TOP_samples <- TOP_samples + load_TOP_samples(TOP_samples_files[i], thin, n.samples)
    }
    # take the average of the samples from multiple partitions
    TOP_samples <- TOP_samples / N
  }

  return(TOP_samples)

}

#' Extract alpha and beta coefficients from TOP posterior samples
#'
#' @param TOP_samples TOP samples combined from all partitions
#' @param tf_cell_combos a table of TF x cell type combinations
#' @param tf_name TF name of interest
#' @param cell_type Cell type of interest
#' @param n.bins Number of DNase or ATAC bins in TOP model (default = 5)
#' @param level The level in the TOP model (bottom, middle, or top),
#' 'bottom' level: TF- and cell-type- specific,
#' 'middle' level: TF-specific, cell-type generic,
#' 'top' level: TF-generic
#'
#' @return
#' @export
#'
extract_TOP_coef_samples <- function(TOP_samples,
                                     tf_cell_combos,
                                     tf_name,
                                     cell_type,
                                     n.bins = 5,
                                     level = c('bottom', 'middle', 'top')){

  level <- match.arg(level)

  TOP_samples <- as.data.frame(TOP_samples)

  if(level == 'bottom' && !missing(tf_name) && !missing(cell_type)){
    tf_cell_combo <- tf_cell_combos[tf_cell_combos$tf_name == tf_name & tf_cell_combos$cell_type == cell_type,]
    tf_id <- tf_cell_combo$tf_id[1]
    cell_id <- tf_cell_combo$cell_id[1]

    # alpha's and beta's are TF- and cell-type-specific model coefficients
    # specific for each TF and cell type combination
    data_id <- sprintf('[%s,%s]', tf_id, cell_id)
    coef_names <- c(paste0('alpha', data_id), paste0('beta', 1:(1+n.bins), data_id))
    coef_samples <- TOP_samples[, coef_names]

  }else if( level == 'middle' && !missing(tf_name) ){
    tf_id <- tf_cell_combos$tf_id[tf_cell_combos$tf_name == tf_name][1]

    # Alpha and Beta's are TF specific, cell-type-generic model coefficients
    data_id <- sprintf('[%s]', tf_id)
    coef_names <- c(paste0('Alpha', data_id), paste0('Beta', 1:(1+n.bins), data_id))
    coef_samples <- TOP_samples[, coef_names]

  }else if(level == 'top'){
    # A and B's are TF-generic model coefficients
    coef_names <- c('A', paste0('B', 1:(1+n.bins)))
    coef_samples <- TOP_samples[, coef_names]

  }

  return(coef_samples)

}

#' Compute TOP posterior mean coefficients for all three levels of TOP model
#'
#' @param TOP_samples TOP samples combined from all partitions
#' @param tf_cell_combos a table of TF x cell type combinations
#' @param n.bins Number of DNase or ATAC bins in TOP model (default = 5)
#'
#' @return
#' @export
#'
extract_TOP_mean_coef <- function(TOP_samples, tf_cell_combos, n.bins = 5){

  # extract bottom level posterior mean coefficients for all TF x cell type combos
  bottom_level_mean_coef <- matrix(NA, nrow = nrow(tf_cell_combos), ncol = 2+n.bins)
  for( i in 1:nrow(tf_cell_combos)){
    bottom_level_coef_samples <- extract_TOP_coef_samples(TOP_samples,
                                                          tf_cell_combos,
                                                          tf_name = tf_cell_combos[i, 'tf_name'],
                                                          cell_type = tf_cell_combos[i, 'cell_type'],
                                                          n.bins = n.bins,
                                                          level = 'bottom')
    bottom_level_mean_coef[i, ] <- colMeans(bottom_level_coef_samples)
  }
  rownames(bottom_level_mean_coef) <- paste(tf_cell_combos$tf_name, tf_cell_combos$cell_type, sep = '.')
  colnames(bottom_level_mean_coef) <- c('Intercept', 'PWM', paste0('Bin', 1:n.bins))

  # extract middle level posterior mean coefficients for all TFs
  tf_name_list <- unique(as.character(tf_cell_combos$tf_name))
  middle_level_mean_coef <- matrix(NA, nrow = length(tf_name_list), ncol = 2+n.bins)
  for( i in 1:length(tf_name_list)){
    middle_level_coef_samples <- extract_TOP_coef_samples(TOP_samples,
                                                          tf_cell_combos,
                                                          tf_name = tf_name_list[i],
                                                          n.bins = n.bins,
                                                          level = 'middle')
    middle_level_mean_coef[i, ] <- colMeans(middle_level_coef_samples)
  }
  rownames(middle_level_mean_coef) <- tf_name_list
  colnames(middle_level_mean_coef) <- c('Intercept', 'PWM', paste0('Bin', 1:n.bins))

  # extract top level posterior mean coefficients
  top_level_coef_samples <- extract_TOP_coef_samples(TOP_samples,
                                                     tf_cell_combos,
                                                     n.bins = n.bins,
                                                     level = 'top')
  top_level_mean_coef <- colMeans(top_level_coef_samples)
  names(top_level_mean_coef) <- c('Intercept', 'PWM', paste0('Bin', 1:n.bins))

  # combine posterior mean coefficients for all three levels
  TOP_mean_coef <- list(top = top_level_mean_coef,
                        middle = middle_level_mean_coef,
                        bottom = bottom_level_mean_coef)

  return(TOP_mean_coef)
}
