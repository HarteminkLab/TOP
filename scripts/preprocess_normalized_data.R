
preprocess_normalized_data <- function(tfs, cell_types, path, pattern,
                                       total_partition = 1, n_part = 1) {
  n_data_read = 0
  # data_files = c()

  if(missing(path)) {
    path = '/home/home5/zhong/xtmp/hierarchical_model_chip/updated_all_data/duke/'
    # old data
    #path = '/home/home5/zhong/xtmp/hierarchical_model_chip/all_data'
  }

  if(missing(pattern)) {
    pattern = '%s_.*_%s_.*M5.*txt$'
    #old pattern for shuffled data
    #pattern = '%s_.*_%s_.*M5.*txt.random.shuffled$'
  }


  combined_data <- list()
  for(tf in tfs) {
    for(cell_type in cell_types) {
      if(total_partition == 1) {
        data_file <- list.files(path = path,
          pattern=pattern %^% c(tf, cell_type), full.names = TRUE)
      } else {
        data_file <- list.files(path = path,
          pattern=pattern %^% c(tf, cell_type), full.names = TRUE)
      }

      if(length(data_file) == 1) {
          # data_files <- c(data_files, data_file)
          data <- read.table(data_file, header = TRUE)
          # I only need the fifth to the 11th columns
          #data <- data[, 5:11]

          # I need columns 5, 7-12
          data <- data[, c(5, 7:12)]
          # rename some columns
          names(data) <- c('pwm', 'dnase.left2_sum', 'dnase.left1_sum',
                           'dnase.motif_sum', 'dnase.right1_sum', 'dnase.right2_sum', 'chip')

          part_length <- as.integer(nrow(data) / total_partition) + 1
          selection <- (1 + (n_part - 1) * part_length):min(nrow(data), n_part*part_length)

          data <- data[selection, ]
          data$tf <- tf
          data$cell_type <- cell_type
          combined_data[[ '%s|%s' %^% c(tf, cell_type) ]] <- data
          n_data_read = n_data_read + 1
      }
      }
    }
  combined_data <- do.call(rbind, combined_data)
  row.names(combined_data) <- NULL
  return(combined_data)
  cat('read', n_data_read, 'datasets\n')
  # return(data_files)
}
