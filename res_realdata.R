library(aricode)

##### Load data set
bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.data.frame(bladder_cancer$X)
labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"

# COSA/PAM Results
i_values <- sprintf("%.1f", seq(0.1, 5, by = 0.1))

for (i in i_values){
  filename <- paste0('data_res/cosapam_bladder_', i, '.rds')
  res <- readRDS(filename)
  cat(i, '\t ARI:', aricode::ARI(res, labels_aggr), '\n')
}

# RPEClust Results
rpeclust_bladder <- readRDS('data_res/rpeclust_bladder.rds')

aricode::ARI(rpeclust_bladder, labels_aggr)

# VarSelLCM Results
varselclust_bladder <- readRDS('data_res/varselcluster_bladder.rds')

aricode::ARI(varselclust_bladder$clustering, labels_aggr)
length(varselclust_bladder$relevant_features)

# PCA/K-Means results
for (i in c(1:2)){
  filename <- paste0('data_res/pca', i, '_kmeans_bladder.rds')
  res <- readRDS(filename)
  cat('Projecting onto the first', i, 'PC(s) gives ARI:', aricode::ARI(res, labels_aggr), '\n')
}

# Sparse PCA/K-Means results
epsilon_vals <- c("1e-03", "1e-02", "1e-01", "1e+00", "2e+00", "5e+00", "1e+01", "2e+01", "5e+01", "1e+02")
for (eps in epsilon_vals){
  filename <- paste0('data_res/spca_bladder_param_', eps, '.rds')
  res <- readRDS(filename)
  cat('1 component, penalty:', eps, '\t ARI:',
      aricode::ARI(res$spca1$clustering, labels_aggr), '\t Non-zero weights:',
      res$spca1$n_nonzero, '\n')
  cat('2 components, penalty:', eps, '\t ARI:',
      aricode::ARI(res$spca2$clustering, labels_aggr), '\t Non-zero weights:',
      res$spca2$n_nonzero_combined, '\n')
}

# Sparse K-Means results
skmeans_bladder_clean_aggr <- readRDS("data_res/skmeans_bladder_clean_aggr.RDS")
aricode::ARI(skmeans_bladder_clean_aggr[[1]]$Cs, labels_aggr)
length(skmeans_bladder_clean_aggr[[1]]$ws)