
#' @title Assemble TOP training data for all TF x cell type combinations,
#' then split training data into 10 partitions
#' @description Prepare the training data for fitting TOP models.
#' It split training data into 10 partitions and
#' assemble training data for all TF x cell type combinations
#' for each of the partitions
#' using the function \code{assemble_partition_training_data}.
#'
#' @param tf_cell_table a data frame listing all TF x cell type combinations
#' and the training data for each combination.
#' It should have at least three columns, with:
#' TF names, cell types, and file names of the individual training data for
#' each TF x cell type combination. The individual training data
#' should be in .rds or text (.txt, or .csv) format.
#' @param logistic_model Logical. If \code{logistic_model = TRUE},
#' prepare assembled data for the logistic version of TOP model.
#' If \code{logistic_model = FALSE}, prepare assembled data for the
#' quantitative occupancy model (default).
#' @param chip_col The column name of ChIP data in the individual training data
#' (default: \dQuote{chip}).
#' @param training_chrs Chromosomes used for training the model
#' (default: odd chromosomes, chr1, chr3, ..., chr21)
#' @param n_partitions Number of partitions to split the training data (default: 10).
#' @param n_cores Number of cores to run in parallel
#' (default: equal to \code{n_partitions}).
#' @param max_sites Max number of candidate sites to keep for
#' each TF x cell type combination (default: 50000). To reduce computation time,
#' randomly select \code{max_sites} candidate sites for
#' each TF x cell type combination, if the number of candidate sites
#' exceeds \code{max_sites}.
#' @param seed Set seed when sampling sites.
#' @return A list of data frames (default: 10),
#' each containing one partition of the
#' training data with all TF x cell type combinations.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores
#'
#' @export
#' @examples
#'
#' # 'tf_cell_table' should have three columns with:
#' # TF names, cell types, and paths to the training data files, like:
#' #  |   tf_name    |   cell_type   |        data_file         |
#' #  |:------------:|:-------------:|:------------------------:|
#' #  |     CTCF     |     K562      |   CTCF.K562.data.rds     |
#' #  |     CTCF     |     A549      |   CTCF.A549.data.rds     |
#' #  |     CTCF     |    GM12878    |   CTCF.GM12878.data.rds  |
#' #  |     ...      |     ...       |   ...                    |
#'
#' # Assemble training data for the quantitative occupancy model,
#' # use odd chromosomes for training, keep at most 50000 candidate sites for
#' # each TF x cell type combination, and split training data into 10 partitions.
#' assembled_training_data <- assemble_training_data(tf_cell_table,
#'                                                   logistic_model = FALSE,
#'                                                   chip_col = 'chip',
#'                                                   training_chrs = paste0('chr', seq(1,21,2)),
#'                                                   n_partitions=10,
#'                                                   max_sites = 50000)
#'
#' # Assemble training data for the logistic version of the model
#' assembled_training_data <- assemble_training_data(tf_cell_table,
#'                                                   logistic_model = TRUE,
#'                                                   chip_col = 'chip_label',
#'                                                   training_chrs = paste0('chr', seq(1,21,2)),
#'                                                   n_partitions=10,
#'                                                   max_sites = 50000)
#'
assemble_training_data <- function(tf_cell_table,
                                   logistic_model=FALSE,
                                   chip_col='chip',
                                   training_chrs=paste0('chr', seq(1,21,2)),
                                   n_partitions=10,
                                   n_cores=n_partitions,
                                   max_sites=50000,
                                   seed=1){

  tf_cell_table <- as.data.frame(tf_cell_table)

  if(ncol(tf_cell_table) < 3){
    stop('tf_cell_table should have at least three columns! ')
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

  # Assemble training data for each partition
  if(missing(n_cores)){
    n_available_cores <- detectCores(logical = FALSE) - 1
    n_cores <- min(n_available_cores, n_partitions)
  }else{
    n_cores <- min(n_cores, n_partitions)
  }

  registerDoParallel(cores=n_cores)

  all_training_data <- foreach(k=1:n_partitions) %dopar% {
    training_data <- assemble_partition_training_data(tf_cell_table,
                                                      logistic_model,
                                                      chip_col,
                                                      training_chrs,
                                                      n_partitions,
                                                      k,
                                                      max_sites/n_partitions,
                                                      seed)
    # Add TF and cell type indices
    training_data$tf_id <- as.integer(factor(training_data$tf_name, levels = tf_list))
    training_data$cell_id <- as.integer(factor(training_data$cell_type, levels = celltype_list))

    training_data
  }

  tf_cell_combos <- unique(all_training_data[[1]][, c('tf_id', 'cell_id', 'tf_name', 'cell_type')])
  cat(nrow(tf_cell_combos), 'TF x cell type combinations assembled in training data. \n')

  return(all_training_data)
}

#' @title Assemble TOP training data for all TF x cell type combinations
#' in a selected partition
#' @description Prepare the training data for fitting TOP models. It first loads
#' the training data for each of the training TF x cell type combinations from
#' \code{tf_cell_table}, then select the candidate sites in training chromosomes,
#' and split the selected candidate sites into partitions (default: 10), and
#' assemble training data for all TF x cell type combinations for a selected partition.
#' @param tf_cell_table A data frame containing the information about
#' the training TFs, cell types and paths of the training data.
#' The first three columns should contain TF name, cell type,
#' and file name containing the corresponding training data.
#' @param logistic_model Logical; if TRUE, use the logistic version of TOP model.
#' @param chip_col The column name of ChIP data in the combined data
#' (default: \dQuote{chip}).
#' @param training_chrs Chromosomes used for training the model
#' Default: odd chromosomes.
#' @param n_partitions Number of partitions to split the training data
#' (default: 10).
#' @param part An integer between 1 and \code{n_partitions} to specify the
#' partition to assemble the training data for.
#' @param max_sites Max number of candidate sites to keep in each partition
#' (default: 5000). Randomly sample \code{max_sites} sites in each partition
#' if the number of candidate sites in each partition exceeds \code{max_sites}.
#' @param seed seed used when sampling sites.
#' @return A data frame of the training data
#' with all TF x cell type combinations in the selected partition.
#'
assemble_partition_training_data <- function(tf_cell_table,
                                             logistic_model=FALSE,
                                             chip_col='chip',
                                             training_chrs=paste0('chr', seq(1,21,2)),
                                             n_partitions=10,
                                             part=1,
                                             max_sites=5000,
                                             seed=1) {

  if(part > n_partitions || part < 1){
    stop('part should be an integer less than n_partitions!')
  }
  cat('Assemble training data for partition', part, '... \n')

  set.seed(seed)

  # Load training data for each TF x cell type combo,
  # split training data into partitions,
  # and select the subset (part) and combine all TF x cell type combinations
  assembled_trainng_data <- list()

  colnames(tf_cell_table)[1:3] <- c('tf_name', 'cell_type', 'data_file')

  for(i in 1:nrow(tf_cell_table)) {

    tf_name <- as.character(tf_cell_table$tf_name[i])
    cell_type <- as.character(tf_cell_table$cell_type[i])
    data_file <- as.character(tf_cell_table$data_file[i])

    if(!file.exists(data_file)){
      if(part == 1){
        message(sprintf('Warning: data of %s in %s cell is not available!\n', tf_name, cell_type))
        cat('Check data file:', data_file, '\n')
      }
      next
    }

    if(grepl('\\.rds$', data_file, ignore.case = TRUE)){
      data <- as.data.frame(readRDS(data_file))
    }else{
      data <- as.data.frame(data.table::fread(data_file))
    }

    if(!chip_col %in% colnames(data)){
      message(sprintf('Warning: data of %s in %s cell does not have %s column!\n',
                      tf_name, cell_type, chip_col))
      next
    }

    if(!'pwm' %in% colnames(data)){
      pwm_col <- grep('pwm', colnames(data), ignore.case = TRUE)
      if(length(pwm_col)==1){
        colnames(data)[pwm_col] <- 'pwm'
      }else{
        stop(sprintf('No column of PWM scores found in data file %s \n', data_file))
      }
    }

    if(length(grep('bin', colnames(data))) == 0){
      bin_cols <- grep('dnase|atac', colnames(data), ignore.case = TRUE)
      if(length(bin_cols) > 0){
        colnames(data)[bin_cols] <- paste0('bin', 1:length(bin_cols))
      }else{
        stop(sprintf('No columns of DNase or ATAC bins found in data file %s \n', data_file))
      }
    }

    if(logistic_model){
      data <- data.frame(data[, -grep('chip', colnames(data))], chip_label = data[, chip_col])
    }else{
      data <- data.frame(data[, -grep('chip', colnames(data))], chip = data[, chip_col])
    }

    # Select the training set
    training_data <- data[data$chr %in% training_chrs,]
    partition_size <- ceiling(nrow(training_data) / n_partitions)
    partition_groups <- as.factor(ceiling(c(1:nrow(training_data))/partition_size))
    data_partitions <- split(training_data, partition_groups)
    training_data_partition <- data.frame(tf_name = tf_name,
                                          cell_type = cell_type,
                                          data_partitions[[part]])

    # Keep a smaller subset of candidate sites in training data
    if(nrow(training_data_partition) > max_sites){
      rows <- sample(1:nrow(training_data_partition), max_sites)
      training_data_partition <- training_data_partition[rows, ]
    }
    assembled_trainng_data[[ paste(tf_name, cell_type, sep = '.') ]] <- training_data_partition

  }

  ## row combine all data
  assembled_trainng_data <- do.call(rbind, assembled_trainng_data)
  row.names(assembled_trainng_data) <- NULL

  return(assembled_trainng_data)

}
