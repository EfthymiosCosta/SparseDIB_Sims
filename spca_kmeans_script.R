library(aricode)
library(parallel)
library(elasticnet)

# Create output directory
if (!dir.exists("spca_res")) dir.create("spca_res")

# Get list of all dataset files
dataset_files <- list.files("datasets", pattern = "^sim_\\d+\\.rds$", full.names = TRUE)
cat("Found", length(dataset_files), "datasets to process\n")

process_dataset <- function(file_path) {
  # Extract sim_id from filename
  sim_id <- as.numeric(gsub(".*sim_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Load dataset
  data <- readRDS(file_path)
  
  # Center the data
  X_centered <- scale(data$X, center = TRUE, scale = FALSE)
  
  # Parameter values to test
  param_vals <- c(1e-3, 1e-2, 1e-1, 1, 2.5, 5, 10, 25, 50, 100)
  
  # Initialize storage for results
  non_zero_vars <- numeric(length(param_vals))
  aris <- numeric(length(param_vals))
  amis <- numeric(length(param_vals))
  
  for (i in seq_along(param_vals)) {
    param_val <- param_vals[i]
    
    # Perform Sparse PCA
    enet_spca <- elasticnet::spca(X_centered,
                                  K = 1, para = param_val,
                                  max.iter = 1000, sparse = "penalty")
    
    # Count non-zero variables
    non_zero_vars[i] <- sum(enet_spca$loadings != 0)
    
    # Project onto first sparse PC
    proj <- X_centered %*% enet_spca$loadings
    
    # Run K-means with true number of clusters
    K <- data$params$n_clust
    kmeans_result <- kmeans(proj, centers = K, nstart = 100)
    
    # Calculate AMI and ARI
    amis[i] <- AMI(kmeans_result$cluster, data$labels)
    aris[i] <- ARI(kmeans_result$cluster, data$labels)
  }
  
  # Save results for this dataset
  result <- list(
    sim_id = sim_id,
    params = data$params,
    param_vals = param_vals,
    non_zero_vars = non_zero_vars,
    aris = aris,
    amis = amis
  )
  
  output_file <- sprintf("spca_res/spca_kmeans_%04d.rds", sim_id)
  saveRDS(result, output_file)
  
  return(sim_id)
}

# Process all datasets in parallel
n_cores <- 100
cat("Processing", length(dataset_files), "datasets using", n_cores, "cores...\n")

cl <- makeCluster(n_cores)
clusterExport(cl, c("process_dataset"))
clusterEvalQ(cl, {
  library(aricode)
  library(elasticnet)
})

results <- parLapply(cl, dataset_files, process_dataset)
stopCluster(cl)

cat("\nSparse PCA K-means analysis complete!\n")
cat("Results saved to: spca_res/\n")