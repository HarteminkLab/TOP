
args <- commandArgs(trailingOnly = TRUE)
total_partition <- as.integer(args[1])
n_part <- as.integer(args[2])
data_type <- args[3]

if(data_type == 'duke') {
  all_tf_cell <- read.table('~/xtmp/hierarchical_model_chip/updated_all_data/stats/table_combined_multi_tfs_OpenChromDnase.txt', header = FALSE)
} else if(data_type == 'uw') {
  all_tf_cell <- read.table('~/xtmp/hierarchical_model_chip/updated_all_data/stats/table_combined_multi_tfs_UwDnase.txt', header = FALSE)
}

tfs <- as.character(unique(all_tf_cell$V1))
cell_types <- as.character(unique(all_tf_cell$V4))

print('reading data')
if(data_type == 'duke') {
  alldata <- preprocess_normalized_data(tfs, cell_types, total_partition=total_partition, n_part=n_part)
} else if(data_type == 'uw') {
  alldata <- preprocess_normalized_data(tfs, cell_types, path = '/home/home5/zhong/xtmp/hierarchical_model_chip/updated_all_data/uw/', total_partition=total_partition, n_part=n_part)
}


alldata$tf <- as.numeric(as.factor(alldata$tf))
alldata$cell_type <- as.numeric(as.factor(alldata$cell_type))

tf_cell <- unique(alldata[, c('tf', 'cell_type')])

alldata$chip <- asinh(alldata$chip)

jags.data <- list(pwm = alldata$pwm, dnase.left2_sum=alldata$dnase.left2_sum,
                  dnase.left1_sum=alldata$dnase.left1_sum,
                  dnase.motif_sum=alldata$dnase.motif_sum,
                  dnase.right1_sum=alldata$dnase.right1_sum,
                  dnase.right2_sum=alldata$dnase.right2_sum,
                  chip=alldata$chip,
                  cell_type=alldata$cell_type,
                  tf=alldata$tf, N=nrow(alldata), n_tfs=length(unique(alldata$tf)),
                  n_cell_types=length(unique(alldata$cell_type))
                  )

save_param = c('A', 'Alpha','alpha',
               'Beta1', 'Beta2', 'Beta3', 'Beta4','Beta5', 'Beta6',
               'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
               'beta1', 'beta2', 'beta3', 'beta4','beta5', 'beta6',
               'tau', 'T', 'TAU')

# result directory
if(data_type == 'duke') {
  directory <- '/home/home5/zhong/xtmp/hierarchical_model_chip/gibbs_sampling_results_updated_data_duke_with_RELA'
} else if(data_type == 'uw') {
  directory <- '/home/home5/zhong/xtmp/hierarchical_model_chip/gibbs_sampling_results_updated_data_uw_with_RELA'
}

print(Sys.info()[['nodename']])
print(Sys.time())
print('init model')
jagsfit <- jags.model(data = jags.data, file='../model/model_all_normalized_data.r', n.adapt = 10000, n.chains = 1);
print(Sys.time())

for(n in 1:100) {
  print(Sys.time())
  cat('sampling %s 10000 samples\n' %^% n)
  jagsfit_samples <- coda.samples(jagsfit, n.iter = 10000, variable.names = save_param, thin = 10);
  file_name <- '%s/ALL_TF_ALL_DATA_totalPartition%s_part%s_adapt10000_%s_%s_samples.rds' %^%
          c(directory, total_partition, n_part, 1 + (n - 1)*10000, n*10000)
  saveRDS(jagsfit_samples, file=file_name)

  print(Sys.time())
}
