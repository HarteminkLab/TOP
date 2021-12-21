
#' @title Assemble TOP training data for all TF x cell type combos,
#' then split training data into partitions
#'
#' @param tf_cell_table a data frame with TF name, cell type,
#' and links to the training data files.
#' @param transform Transform of ChIP read counts.
#' Options are 'asinh' (default), 'log2', 'sqrt', 'none'.
#' @param training_chrs Chromosomes used for training the model.
#' @param sites_limit Max number of candidate sites in each partition.
#' @param total_partitions split the data into partitions (default = 10)
#' to run Gibbs sampling in parallel.
#' @param k which partition to use for training the model
#'
#' @return Returns a data frame of training data with all TFs and cell type combos.
#' @export
#'
assemble_TOP_training_data <- function(tf_cell_table,
                                       transform = c('asinh', 'log2', 'sqrt', 'none'),
                                       training_chrs = paste0('chr', seq(1,21,2)),
                                       sites_limit = Inf,
                                       total_partitions = 10,
                                       k = 1) {

  cat('Split and combine training data for partition', k, '... \n')

  transform <- match.arg(transform)
  cat(transform, 'transform ChIP counts. \n')

  # Load training data for each TF x cell type combo,
  # split training data into partitions,
  # and select the subset (k) and combine all TF x cell type combos
  assembled_trainng_data <- list()
  summary_training_table <- NULL
  data_count <- 0

  for(i in 1:nrow(tf_cell_table)) {

    tf_name <- as.character(tf_cell_table$tf_name[i])
    cell_type <- as.character(tf_cell_table$cell_type[i])
    data_file <- as.character(tf_cell_table$data_file[i])

    if(!file.exists(data_file) || file.size(data_file) == 0){
      if(k == 1){
        cat('Warning:', tf_name, 'in', cell_type, 'data file does not exist or is empty!\n')
        cat('Check:', data_file, '\n')
      }
      next
    }

    data <- as.data.frame(readRDS(data_file))

    # Transform the chip data
    if (transform == 'asinh') {
      data$chip <- asinh(data$chip)
    } else if (transform == 'log2') {
      data$chip <- log2(data$chip + 1)
    } else if (transform == 'sqrt') {
      data$chip <- sqrt(data$chip)
    }

    # Select training and test set
    training_data <- data[data$chr %in% training_chrs,]

    partition_size <- ceiling(nrow(training_data) / total_partitions)
    partition_groups <- as.factor(ceiling(c(1:nrow(training_data))/partition_size))
    data_partitions <- split(training_data, partition_groups)

    training_data_partition <- data.frame(tf_name = tf_name, cell_type = cell_type, data_partitions[[k]])
    # set the limit of training candidate sites
    if(nrow(training_data_partition) > sites_limit){
      training_data_partition <- training_data_partition[sample(1:nrow(training_data_partition), sites_limit), ]
    }

    assembled_trainng_data[[ paste(tf_name, cell_type, sep = '|') ]] <- training_data_partition

    summary_training_table <- rbind(summary_training_table,
                                    data.frame(tf_name = tf_name, cell_type = cell_type, n_training_sites = nrow(training_data_partition)))

    data_count <- data_count + 1

  }

  ## row combine all data
  assembled_trainng_data <- do.call(rbind, assembled_trainng_data)
  row.names(assembled_trainng_data) <- NULL

  cat('Assembled', data_count, 'datasets. \n')

  return(assembled_trainng_data)

}


#' @title Assemble TOP (logistic version) training data for all TF x cell type combos,
#' then split training data into partitions
#' #'
#' @param tf_cell_table a data frame with TF name, cell type,
#' and links to the training data files.
#' @param chiplabel_colname The column name of ChIP peak label in the combined data.
#' @param training_chrs Chromosomes used for training the model.
#' @param sites_limit Max number of candidate sites in each partition.
#' @param total_partitions split the data into partitions (default = 10)
#' to run Gibbs sampling in parallel.
#' @param k which partition to use for training the model
#'
#' @return Returns a data frame of training data with all TFs and cell type combos.
#' @export
#'
assemble_TOP_logistic_training_data <- function(tf_cell_table,
                                                chiplabel_colname = 'chip_label',
                                                training_chrs = paste0('chr', seq(1,21,2)),
                                                sites_limit = Inf,
                                                total_partitions = 10,
                                                k = 1) {

  cat('Split and combine training data for partition', k, '... \n')

  # Load training data for each TF x cell type combo,
  # split training data into partitions,
  # and select the subset (k) and combine all TF x cell type combos
  assembled_trainng_data <- list()
  summary_training_table <- NULL
  data_count <- 0

  for(i in 1:nrow(tf_cell_table)) {

    tf_name <- as.character(tf_cell_table$tf_name[i])
    cell_type <- as.character(tf_cell_table$cell_type[i])
    data_file <- as.character(tf_cell_table$data_file[i])

    if(!file.exists(data_file) || file.size(data_file) == 0){
      if(k == 1){
        cat('Warning:', tf_name, 'in', cell_type, 'data file does not exist or is empty!\n')
        cat('Check:', data_file, '\n')
      }
      next
    }

    data <- as.data.frame(readRDS(data_file))

    # cat('ChIP occupancy column:' , chip_colname, '\n')
    if(!chiplabel_colname %in% colnames(data)){
      cat('Warning:', tf_name, 'in', cell_type, 'data file does not have', chiplabel_colname, 'column!\n')
      next
    }

    data <- data.frame(data[, -grep('chip', colnames(data))], chip_label = data[, chiplabel_colname])

    # Select training and test set
    training_data <- data[data$chr %in% training_chrs,]

    partition_size <- ceiling(nrow(training_data) / total_partitions)
    partition_groups <- as.factor(ceiling(c(1:nrow(training_data))/partition_size))
    data_partitions <- split(training_data, partition_groups)

    training_data_partition <- data.frame(tf_name = tf_name, cell_type = cell_type, data_partitions[[k]])
    # set the limit of training candidate sites
    if(nrow(training_data_partition) > sites_limit){
      training_data_partition <- training_data_partition[sample(1:nrow(training_data_partition), sites_limit), ]
    }

    assembled_trainng_data[[ paste(tf_name, cell_type, sep = '|') ]] <- training_data_partition

    summary_training_table <- rbind(summary_training_table,
                                    data.frame(tf_name = tf_name, cell_type = cell_type, n_training_sites = nrow(training_data_partition)))

    data_count <- data_count + 1

  }

  ## row combine all data
  assembled_trainng_data <- do.call(rbind, assembled_trainng_data)
  row.names(assembled_trainng_data) <- NULL

  cat('Assembled', data_count, 'datasets. \n')

  return(assembled_trainng_data)

}
