if (!dir.exists("data_res")) dir.create("data_res")

# Runs sparse K Means determining the value of wbound based on GAP statistic
# Returns the final clustering outpuT
SparseKMeans_full <- function(X, labels, wbounds, nperms, nstart, maxiter){
  library(sparcl)
  perm_test <- KMeansSparseCluster.permute(x = X,
                                           K = length(unique(labels)),
                                           nperms = nperms,
                                           wbounds = wbounds)
  skmeans_res <- KMeansSparseCluster(x = X,
                                     K = length(unique(labels)),
                                     wbounds = perm_test$bestw,
                                     nstart = nstart,
                                     maxiter = maxiter)
  return(skmeans_res)
}

# Bladder cancer dataset
bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"
X <- as.data.frame(bladder_cancer$X)

skmeans_res_bladder_aggr <- SparseKMeans_full(X = X,
                                         labels = labels_aggr,
                                         wbounds = seq(2, 101, length.out = 100),
                                         nperms = 100,
                                         nstart = 100,
                                         maxiter = 100)
saveRDS(skmeans_res_bladder_aggr, file = 'data_res/skmeans_bladder_clean_aggr.RDS')