library(aricode)
library(parallel)
library(elasticnet)

# Create output directory
if (!dir.exists("data_res")) dir.create("data_res")
bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.data.frame(bladder_cancer$X)
labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"
K <- length(unique(labels_aggr))
X_centered <- scale(X, center = TRUE, scale = FALSE)

# Parameter values to test
param_vals <- c(0.001, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)

cat("Running Sparse PCA + K-Means with", length(param_vals), "parameter values...\n")
cat("Using", K, "clusters\n")

process_param <- function(param_val) {
  
  # Sparse PCA with K=1
  enet_spca1 <- elasticnet::spca(X_centered,
                                 K = 1, para = param_val,
                                 max.iter = 1000, sparse = "penalty")
  
  # Non-zero loading indices/names
  nonzero_idx1 <- which(enet_spca1$loadings != 0)
  nonzero_names1 <- colnames(X)[nonzero_idx1]
  
  # Project onto first sparse PC
  proj1 <- X_centered %*% enet_spca1$loadings
  
  # Run K-means
  kmeans_result1 <- kmeans(proj1, centers = K, nstart = 100)
  
  # Sparse PCA with K=2
  enet_spca2 <- elasticnet::spca(X_centered,
                                 K = 2, para = param_val,
                                 max.iter = 1000, sparse = "penalty")
  
  # Non-zero loading indices/names for each PC
  nonzero_idx2_pc1 <- which(enet_spca2$loadings[, 1] != 0)
  nonzero_idx2_pc2 <- which(enet_spca2$loadings[, 2] != 0)
  nonzero_names2_pc1 <- colnames(X)[nonzero_idx2_pc1]
  nonzero_names2_pc2 <- colnames(X)[nonzero_idx2_pc2]
  
  # Combined non-zero variables across both PCs
  nonzero_idx2_combined <- unique(c(nonzero_idx2_pc1, nonzero_idx2_pc2))
  nonzero_names2_combined <- colnames(X)[nonzero_idx2_combined]
  
  # Project onto first 2 sparse PCs
  proj2 <- X_centered %*% enet_spca2$loadings
  
  # Run K-means
  kmeans_result2 <- kmeans(proj2, centers = K, nstart = 100)
  
  result <- list(
    param_val = param_val,
    spca1 = list(
      clustering = kmeans_result1$cluster,
      nonzero_indices = nonzero_idx1,
      nonzero_names = nonzero_names1,
      n_nonzero = length(nonzero_idx1)
    ),
    spca2 = list(
      clustering = kmeans_result2$cluster,
      nonzero_indices_pc1 = nonzero_idx2_pc1,
      nonzero_names_pc1 = nonzero_names2_pc1,
      nonzero_indices_pc2 = nonzero_idx2_pc2,
      nonzero_names_pc2 = nonzero_names2_pc2,
      nonzero_indices_combined = nonzero_idx2_combined,
      nonzero_names_combined = nonzero_names2_combined,
      n_nonzero_pc1 = length(nonzero_idx2_pc1),
      n_nonzero_pc2 = length(nonzero_idx2_pc2),
      n_nonzero_combined = length(nonzero_idx2_combined)
    )
  )
  
  output_file <- sprintf("data_res/spca_bladder_param_%.0e.RDS", param_val)
  saveRDS(result, output_file)
  
  return(param_val)
}

n_cores <- 10
cl <- makeCluster(n_cores)
clusterExport(cl, c("X_centered", "K", "X"))
clusterEvalQ(cl, {
  library(elasticnet)
})

results <- parLapply(cl, param_vals, process_param)
stopCluster(cl)

cat("\nSparse PCA K-means analysis complete!\n")