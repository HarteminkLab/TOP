
# Load TOP posterior samples
load_TOP_samples <- function(sample_file, thin = 1, n_samples = 1000) {

  if (!file.exists(sample_file)) {
    stop(paste("Sample file: ", sample_file, "is not available!"))
  }
  samples <- readRDS(sample_file)

  # merge multiple MCMC chains
  if ( is.list(samples) ){
    samples <- as.data.frame(do.call(rbind, samples))
  }

  # cat('Load', nrow(samples), 'samples... \n')

  if (thin > 1) {
    samples <- samples[seq(from = 1, to = nrow(samples), by = thin), ]
  } else if (nrow(samples) > n_samples) {
    # cat('Randomly choose', n_samples, 'samples... \n')
    samples <- samples[seq(from = 1, to = nrow(samples), length.out = n_samples),]
  }

  return(samples)

}

# Combine (average) TOP samples from multiple partitions
combine_TOP_samples <- function(sample_files, thin = 1, n_samples = 1000) {

  N <- length(sample_files)
  cat('Combining', length(sample_files), 'sample partitions ... \n')
  samples <- load_TOP_samples(sample_files[1], thin, n_samples)

  if(N > 1) {
    for(i in 2:N) {
      samples <- samples + load_TOP_samples(sample_files[i], thin, n_samples)
    }
    # take the average of the samples from multiple partitions
    samples <- samples / N
  }

  return(samples)

}

## extract coefficients from alpha's and beta's
extract_TOP_coef_samples <- function(TOP_samples,
                                     tf_cell_combos,
                                     tf_name,
                                     cell_type,
                                     n_bins = 5,
                                     level = c('top', 'middle', 'bottom')){

  level <- match.arg(level)
  TOP_samples <- as.data.frame(TOP_samples)

  if(level == 'bottom' && !missing(tf_name) && !missing(cell_type)){
    tf_cell_combo <- tf_cell_combos[tf_cell_combos$tf_name == tf_name & tf_cell_combos$cell_type == cell_type,]
    tf_id <- tf_cell_combo$tf_id[1]
    cell_id <- tf_cell_combo$cell_id[1]

    ## alpha's and beta's are TF and cell type specific model coefficients for each TF and cell type combination
    data_id <- sprintf('[%s,%s]', tf_id, cell_id)
    coef_names <- c(paste0('alpha', data_id), paste0('beta', 1:(1+n_bins), data_id))
    coef_samples <- TOP_samples[, coef_names]

  }else if( level == 'middle' && !missing(tf_name) ){
    tf_id <- tf_cell_combos$tf_id[tf_cell_combos$tf_name == tf_name][1]

    ## Alpha and Beta's are TF specific (cell type generic) model coefficients for each TF
    data_id <- sprintf('[%s]', tf_id)
    coef_names <- c(paste0('Alpha', data_id), paste0('Beta', 1:(1+n_bins), data_id))
    coef_samples <- TOP_samples[, coef_names]

  }else if(level == 'top'){
    ## A and B's are TF generic model coefficients
    coef_names <- c('A', paste0('B', 1:(1+n_bins)))
    coef_samples <- TOP_samples[, coef_names]

  }

  return(coef_samples)

}
