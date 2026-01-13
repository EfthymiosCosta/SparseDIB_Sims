library(IBclust)
library(parallel)
library(aricode)

# Source required functions
source('src/sparse_DIB.R')
source('src/mi_comps.R')
source('src/dykstra.R')
source('src/coord_to_pxy_weights.R')

# Create output directory if it doesn't exist
if (!dir.exists('sdib_res')) {
  dir.create('sdib_res', recursive = TRUE)
}

sparsity_values <- seq(0.4, 10, by = 0.2)

# Function to process a single dataset
process_dataset <- function(sim_id) {
  # Format sim_id with leading zeros
  sim_id_str <- sprintf("%04d", sim_id)
  
  # Load data
  data <- readRDS(paste0('datasets/sim_', sim_id_str, '.rds'))
  X <- as.data.frame(data$X)
  labels <- data$labels
  
  # Initialize result storage
  ari_vec <- c()
  ami_vec <- c()
  entropy_vec <- c()
  nonzero_weights_vec <- c()
  
  # Loop over sparsity values
  for (sparsity_val in sparsity_values) {
    res_temp <- sparse_DIB(X, 
                           ncl = data$params$n_clust, 
                           s = -1,
                           contkernel = "gaussian",
                           nomkernel = "aitchisonaitken",
                           ordkernel = "liracine",
                           cat_first = FALSE,
                           nstart = 100,
                           maxiter = 100,
                           maxsteps = 15,
                           randinit = NULL,
                           verbose = FALSE,
                           sparsity = sparsity_val,
                           mRMR = FALSE,
                           mi_matrix = NULL,
                           n_cores = 1) 
    
    ari <- ARI(res_temp$best_clust$Cluster, labels)
    ami <- AMI(res_temp$best_clust$Cluster, labels)
    normalized_weights <- res_temp$weights / sum(res_temp$weights)
    entr <- IBclust:::entropy(normalized_weights)
    n_nonzero <- sum(res_temp$weights > 0)
    
    ari_vec <- c(ari_vec, ari)
    ami_vec <- c(ami_vec, ami)
    entropy_vec <- c(entropy_vec, entr)
    nonzero_weights_vec <- c(nonzero_weights_vec, n_nonzero)
  }
  
  results <- list(
    sim_id = sim_id,
    params = data$params,
    sparsity_values = sparsity_values,
    ari = ari_vec,
    ami = ami_vec,
    entropy = entropy_vec,
    n_nonzero_weights = nonzero_weights_vec
  )
  
  saveRDS(results, paste0('sdib_res/sdib_', sim_id_str, '.rds'))
  
  return(results)
}

all_sim_ids <- 1:9600

# Parallel processing across datasets
cl <- makeCluster(123)

clusterExport(cl, c("process_dataset", "sparse_DIB", "coord_to_pxy_weights", 
                    "compute_py_x_single", "dykstra_projection", "sparsity_values"))
clusterEvalQ(cl, {
  library(IBclust)
  library(aricode)
  source('src/mi_comps.R')
  source('src/dykstra.R')
})

# Run in parallel
results <- parLapply(cl, all_sim_ids, process_dataset)

stopCluster(cl)

cat("All datasets processed successfully!\n")