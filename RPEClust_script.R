library(RPEClust)
library(aricode)
library(parallel)

# Create output directory
if (!dir.exists("rpeclust_res")) dir.create("rpeclust_res")

# Get list of all dataset files
dataset_files <- list.files("datasets", pattern = "^sim_\\d+\\.rds$", full.names = TRUE)
cat("Found", length(dataset_files), "datasets to process\n")

process_dataset <- function(file_path) {
  # Extract sim_id from filename
  sim_id <- as.numeric(gsub(".*sim_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Load dataset
  data <- readRDS(file_path)
  
  # Run RPGMMClu
  K <- data$params$n_clust
  
  res <- RPGMMClu(x = data$X, true.cl = data$labels,
                  g = K, B = 1000)
  
  # Calculate AMI and ARI
  ami_score <- AMI(res$ensemble$label.vec, data$labels)
  ari_score <- ARI(res$ensemble$label.vec, data$labels)
  
  # Save results
  result <- list(
    sim_id = sim_id,
    ami = ami_score,
    ari = ari_score,
    params = data$params
  )
  
  output_file <- sprintf("rpeclust_res/rpeclust_%04d.rds", sim_id)
  saveRDS(result, output_file)
  
  return(sim_id)
}

# Setup parallel processing
n_cores <- 123
cl <- makeCluster(n_cores)

# Export necessary libraries and functions to cluster
clusterEvalQ(cl, {
  library(RPEClust)
  library(aricode)
})

# Export the process_dataset function to cluster
clusterExport(cl, "process_dataset")

cat("Processing datasets in parallel with", n_cores, "cores...\n")

# Process all datasets in parallel
results <- parLapply(cl, dataset_files, process_dataset)

# Stop cluster
stopCluster(cl)

cat("\nRPEClust analysis complete!\n")
cat("Results saved to: rpeclust_res/\n")