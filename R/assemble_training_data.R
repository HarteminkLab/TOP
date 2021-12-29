#' @title Assemble TOP training data for all TF x cell type combos in one partition
#'
#' @param tf_cell_table a data frame with at least three columns:
#' TF name, cell type, and file name containing the corresponding training data.
#' @param logistic.model Logical; if TRUE, use the logistic version of TOP model.
#' @param chip_colname The column name of ChIP data in the combined data (default: "chip").
#' @param training_chrs Chromosomes used for training the model (default: odd chromosomes)
#' @param n.partitions Total number of partitions to split the training data (default: 10).
#' @param part specifies which partition to assemble the training data for
#' @param max.sites Max number of candidate sites in each partition (default: 50000/n.partitions).
#' @param seed seed used when sampling sites (default: 123).
#'
#' @export
#'
assemble_partition_training_data <- function(tf_cell_table,
                                             logistic.model = FALSE,
                                             chip_colname = 'chip',
                                             training_chrs = paste0('chr', seq(1,21,2)),
                                             n.partitions = 10,
                                             part,
                                             max.sites = 50000/n.partitions,
                                             seed = 123) {

  cat('Assemble training data for partition', part, '... \n')

  set.seed(seed)

  # Load training data for each TF x cell type combo,
  # split training data into partitions,
  # and select the subset (part) and combine all TF x cell type combos
  assembled_trainng_data <- list()

  colnames(tf_cell_table)[1:3] <- c('tf_name', 'cell_type', 'data_file')

  for(i in 1:nrow(tf_cell_table)) {

    tf_name <- as.character(tf_cell_table$tf_name[i])
    cell_type <- as.character(tf_cell_table$cell_type[i])
    data_file <- as.character(tf_cell_table$data_file[i])

    if(!file.exists(data_file)){
      if(part == 1){
        cat('Warning:', tf_name, 'in', cell_type, 'data file is not available!\n')
        cat('Check data file:', data_file, '\n')
      }
      next
    }

    if(grepl('.rds', data_file, ignore.case = TRUE)){
      data <- as.data.frame(readRDS(data_file))
    }else{
      data <- as.data.frame(data.table::fread(data_file))
    }


    if(!chip_colname %in% colnames(data)){
      cat('Warning:', tf_name, 'in', cell_type, 'data file does not have', chip_colname, 'column!\n')
      next
    }

    if(!'pwm' %in% colnames(data)){
       pwm_col <- grep('pwm', colnames(data), ignore.case = TRUE)
       if(length(pwm_col)==1){
         colnames(data)[pwm_col] <- 'pwm'
       }else{
         stop('No \"pwm\" column found in data file! \n')
       }
    }

    if(length(grep('bin', colnames(data))) == 0){
      bin_cols <- grep('dnase|atac', colnames(data), ignore.case = TRUE)
      if(length(bin_cols) > 0){
        colnames(data)[bin_cols] <- paste0('bin', 1:length(bin_cols))
      }else{
        stop('No \"pwm\" column found in data file! \n')
      }
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
    assembled_trainng_data[[ paste(tf_name, cell_type, sep = '.') ]] <- training_data_partition

  }

  ## row combine all data
  assembled_trainng_data <- do.call(rbind, assembled_trainng_data)
  row.names(assembled_trainng_data) <- NULL

  return(assembled_trainng_data)

}

#' @title Assemble TOP training data for all TF x cell type combos,
#' then split training data into 10 partitions
#'
#' @param tf_cell_table_file a tab delimited file with at least three columns:
#' The first three columns should be TF names, cell types, and the training data files.
#' @param logistic.model Logical; if TRUE, use the logistic version of TOP model.
#' @param chip_colname The column name of ChIP data in the combined data (default: "chip").
#' @param training_chrs Chromosomes used for training the model (default: odd chromosomes)
#' @param n.partitions Total number of partitions to split the training data (default: 10).
#' @param n.cores Number of cores to use in parallel
#' (default: equal to the number of partitions).
#' @param max.sites Max number of candidate sites in each partition (default: 50000/n.partitions).
#' @param seed seed used when sampling sites (default: 123).
#' @import doParallel
#' @import foreach
#'
#' @export
#'
assemble_TOP_training_data <- function(tf_cell_table_file,
                                       logistic.model=FALSE,
                                       chip_colname='chip',
                                       training_chrs=paste0('chr', seq(1,21,2)),
                                       n.partitions=10,
                                       n.cores=n.partitions,
                                       max.sites=50000,
                                       seed=123){

  tf_cell_table <- as.data.frame(data.table::fread(tf_cell_table_file))

  if(ncol(tf_cell_table) < 3){
    stop('The table should have at least three columns separated by tab! ')
  }
  colnames(tf_cell_table)[1:3] <- c('tf_name', 'cell_type', 'data_file')

  cat('Assemble TOP training data...\n')
  tf_list <- sort(unique(as.character(tf_cell_table$tf_name)))
  celltype_list <- sort(unique(as.character(tf_cell_table$cell_type)))
  cat('TFs:', tf_list, '\n')
  cat('Cell types:', celltype_list, '\n')

  tf_cell_table$tf_name <- factor(tf_cell_table$tf_name, levels = tf_list)
  tf_cell_table$cell_type <- factor(tf_cell_table$cell_type, levels = celltype_list)
  tf_cell_table <- tf_cell_table[with(tf_cell_table, order(tf_name, cell_type)),]

  registerDoParallel(cores=n.cores)
  cat('Using', getDoParWorkers(), 'cores in parallel. \n')

  all_training_data <- foreach(k=1:n.partitions) %dopar% {
    training_data <- assemble_partition_training_data(tf_cell_table,
                                                      logistic.model,
                                                      chip_colname,
                                                      training_chrs,
                                                      n.partitions,
                                                      k,
                                                      max.sites/n.partitions,
                                                      seed)
    # Add TF and cell type indices
    training_data$tf_id <- as.integer(factor(training_data$tf_name, levels = tf_list))
    training_data$cell_id <- as.integer(factor(training_data$cell_type, levels = celltype_list))

    training_data
  }

  # Save a table listing all TF and cell type combinations
  tf_cell_combos <- unique(all_training_data[[1]][, c('tf_id', 'cell_id', 'tf_name', 'cell_type')])
  cat(nrow(tf_cell_combos), 'TF x cell type combos assembled in training data. \n')

  return(all_training_data)
}

