#' @title Assemble TOP training data for all TF x cell type combos in one partition
#'
#' @param tf_cell_table a data frame with TF name, cell type,
#' and links to the training data files.
#' @param logistic.model If TRUE, use logistic version of the model
#' @param chip_colname The column name of ChIP data in the combined data.
#' @param training_chrs Chromosomes used for training the model.
#' @param n.partitions split the data into partitions (default = 10)
#' to run Gibbs sampling in parallel.
#' @param part which partition to use for training the model
#' @param max.sites Max number of candidate sites in each partition.
#' @param seed seed used when sampling sites.
#'
#' @return Returns a data frame of training data with all TFs and cell type combos.
#' @export
#'
assemble_partition_training_data <- function(tf_cell_table,
                                             logistic.model = FALSE,
                                             chip_colname = 'chip',
                                             training_chrs = paste0('chr', seq(1,21,2)),
                                             n.partitions = 10,
                                             part,
                                             max.sites = 10000,
                                             seed = 1) {

  cat('Assemble training data for partition', part, '... \n')

  set.seed(seed)

  # Load training data for each TF x cell type combo,
  # split training data into partitions,
  # and select the subset (part) and combine all TF x cell type combos

  assembled_trainng_data <- list()
  data_count <- 0

  for(i in 1:nrow(tf_cell_table)) {

    tf_name <- as.character(tf_cell_table$tf_name[i])
    cell_type <- as.character(tf_cell_table$cell_type[i])
    data_file <- as.character(tf_cell_table$data_file[i])

    if(!file.exists(data_file)){
      if(part == 1){
        cat('Warning:', tf_name, 'in', cell_type, 'data file is not available!\n')
        cat('Check:', data_file, '\n')
      }
      next
    }

    data <- as.data.frame(readRDS(data_file))

    if(!chip_colname %in% colnames(data)){
      cat('Warning:', tf_name, 'in', cell_type, 'data file does not have', chip_colname, 'column!\n')
      next
    }

    if(logistic.model){
      data <- data.frame(data[, -grep('chip', colnames(data))], chip_label = data[, chip_colname])
    }else{
      data <- data.frame(data[, -grep('chip', colnames(data))], chip = data[, chip_colname])
    }

    # Select training and test set
    training_data <- data[data$chr %in% training_chrs,]

    partition_size <- ceiling(nrow(training_data) / n.partitions)
    partition_groups <- as.factor(ceiling(c(1:nrow(training_data))/partition_size))
    data_partitions <- split(training_data, partition_groups)

    training_data_partition <- data.frame(tf_name = tf_name, cell_type = cell_type, data_partitions[[part]])
    # set the limit of training candidate sites
    if(nrow(training_data_partition) > max.sites){
      training_data_partition <- training_data_partition[sample(1:nrow(training_data_partition), max.sites), ]
    }
    assembled_trainng_data[[ paste(tf_name, cell_type, sep = '|') ]] <- training_data_partition

    data_count <- data_count + 1

  }

  ## row combine all data
  assembled_trainng_data <- do.call(rbind, assembled_trainng_data)
  row.names(assembled_trainng_data) <- NULL

  cat('Assembled', data_count, 'datasets. \n')

  return(assembled_trainng_data)

}

#' @title Assemble TOP training data for all TF x cell type combos,
#' then split training data into 10 partitions
#'
#' @param tf_cell_table a data frame with TF name, cell type,
#' and links to the training data files.
#' @param training_data_dir Directory for saving training data
#' @param training_data_name Prefix for training data file names
#' @param logistic.model If TRUE, use logistic version of the model
#' @param chip_colname The column name of ChIP data in the combined data.
#' @param training_chrs Chromosomes used for training the model.
#' @param n.partitions split the data into partitions (default = 10)
#' to run Gibbs sampling in parallel.
#' @param max.sites Max number of candidate sites in each partition.
#' @param seed seed used when sampling sites.
#' @import doParallel
#' @import foreach
#' @importFrom data.table fwrite
#'
#' @export
#'
assemble_TOP_training_data <- function(tf_cell_table,
                                       training_data_dir="./",
                                       training_data_name='TOP_training_data',
                                       logistic.model=FALSE,
                                       chip_colname="chip",
                                       training_chrs=paste0('chr', seq(1,21,2)),
                                       n.partitions=10,
                                       max.sites=50000,
                                       seed=1){

  if(!dir.exists(training_data_dir)){
    dir.create(training_data_dir, showWarnings = FALSE, recursive = TRUE)
  }

  cat("Assemble TOP training data...\n")
  tf_list <- sort(unique(tf_cell_table$tf_name))
  celltype_list <- sort(unique(tf_cell_table$cell_type))
  cat("TFs:", tf_list, "\n")
  cat("Cell types:", celltype_list, "\n")

  tf_cell_table$tf_name <- factor(tf_cell_table$tf_name, levels = tf_list)
  tf_cell_table$cell_type <- factor(tf_cell_table$cell_type, levels = celltype_list)
  tf_cell_table <- tf_cell_table[with(tf_cell_table, order(tf_name, cell_type)),]

  doParallel::registerDoParallel(cores=n.partitions)
  cat("Using", foreach::getDoParWorkers(), "cores in parallel. \n")

  all_training_data <- foreach(k=1:n.partitions) %dopar% {
    training_data <- assemble_partition_training_data(tf_cell_table,
                                                      chip_colname,
                                                      logistic.model,
                                                      training_chrs,
                                                      n.partitions,
                                                      k,
                                                      max.sites/n.partitions,
                                                      seed)
    # Add TF and cell type indices
    training_data$tf_id <- as.integer(factor(training_data$tf_name, levels = tf_list))
    training_data$cell_id <- as.integer(factor(training_data$cell_type, levels = celltype_list))
    saveRDS(training_data,
            file.path(training_data_dir, paste0(training_data_name, '.partition', k, '.rds')))
    training_data
  }

  # Save a table with all TF and cell type combinations
  tf_cell_combos <- unique(training_data[, c('tf_id', 'cell_id', 'tf_name', 'cell_type')])
  cat(nrow(tf_cell_combos), 'TF x cell type combos assembled in training data. \n')
  fwrite(tf_cell_combos,
         file.path(training_data_dir, paste0(training_data_name, '_tf_cell_combos.txt')),
         sep = '\t')

  return(all_training_data)
}

