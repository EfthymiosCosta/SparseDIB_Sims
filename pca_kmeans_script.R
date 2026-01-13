library(aricode)
library(parallel)

# Create output directory
if (!dir.exists("pca_res")) dir.create("pca_res")

# Get list of all dataset files
dataset_files <- list.files("datasets", pattern = "^sim_\\d+\\.rds$", full.names = TRUE)

cat("Found", length(dataset_files), "datasets to process\n")

process_dataset <- function(file_path) {
  # Extract sim_id from filename
  sim_id <- as.numeric(gsub(".*sim_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Load dataset
  data <- readRDS(file_path)
  
  # Perform PCA
  pca_result <- prcomp(data$X, center = TRUE, scale. = TRUE)
  
  # Project onto first PC
  pc1_projection <- pca_result$x[, 1, drop = FALSE]
  
  # Run K-means with true number of clusters
  K <- data$params$n_clust
  kmeans_result <- kmeans(pc1_projection, centers = K, nstart = 100)
  
  # Calculate AMI and ARI
  ami_score <- AMI(kmeans_result$cluster, data$labels)
  ari_score <- ARI(kmeans_result$cluster, data$labels)
  
  # Save results
  result <- list(
    sim_id = sim_id,
    ami = ami_score,
    ari = ari_score,
    params = data$params
  )
  
  output_file <- sprintf("pca_res/pca_kmeans_%04d.rds", sim_id)
  saveRDS(result, output_file)
  
  return(sim_id)
}

# Process all datasets in parallel
n_cores <- 100
cat("Processing", length(dataset_files), "datasets using", n_cores, "cores...\n")

cl <- makeCluster(n_cores)
clusterExport(cl, c("process_dataset"))
clusterEvalQ(cl, library(aricode))

results <- parLapply(cl, dataset_files, process_dataset)

stopCluster(cl)

cat("\nPCA K-means analysis complete!\n")
cat("Results saved to: pca_res/\n")