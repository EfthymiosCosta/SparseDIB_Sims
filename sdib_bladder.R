if (!dir.exists("data_res")) dir.create("data_res")
library(IBclust)
library(parallel)

# Source required functions
source('src/sparse_DIB.R')
source('src/mi_comps.R')
source('src/dykstra.R')
source('src/coord_to_pxy_weights.R')

sparsity_values <- c(seq(0.1, 5.9, by = 0.1), seq(6, 100, by = 1))

# Bladder cancer dataset
bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.data.frame(bladder_cancer$X)
labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"
km_test <- kmeans(X, centers = length(unique(labels_aggr)), nstart = 1000)
cat("K-Means ARI:", aricode::ARI(km_test$cluster, labels_aggr), "\n")

# Define function to run for each sparsity value
run_sparse_DIB <- function(sparsity_val) {
  sdib_res <- sparse_DIB(X, 
                         ncl = length(unique(labels_aggr)), 
                         s = -1,
                         contkernel = "gaussian",
                         nomkernel = "aitchisonaitken",
                         ordkernel = "liracine",
                         cat_first = FALSE,
                         nstart = 100,
                         maxiter = 100,
                         maxsteps = 50,
                         randinit = NULL,
                         verbose = FALSE,
                         sparsity = sparsity_val,
                         init_weights = "warm",
                         n_cores = 1)
  
  filename <- paste0("data_res/sdib_km_res_bladder_clean_aggr_", sparsity_val, ".RDS")
  saveRDS(sdib_res, file = filename)
}

# Create cluster
cl <- makeCluster(50)

# Export necessary objects and functions to each worker
clusterExport(cl, c("X", "labels_aggr", "sparse_DIB"))

# Load required libraries and source files on each worker
clusterEvalQ(cl, {
  library(IBclust)
  source('src/sparse_DIB.R')
  source('src/mi_comps.R')
  source('src/dykstra.R')
  source('src/coord_to_pxy_weights.R')
})

# Run in parallel
parLapply(cl, sparsity_values, run_sparse_DIB)

# Stop cluster
stopCluster(cl)