library(VarSelLCM)
library(aricode)
library(parallel)

# Create output directory
if (!dir.exists("varselcluster_res")) dir.create("varselcluster_res")

# Get list of all dataset files
dataset_files <- list.files("datasets", pattern = "^sim_\\d+\\.rds$", full.names = TRUE)

cat("Found", length(dataset_files), "datasets to process\n")

process_dataset <- function(file_path) {
  # Extract sim_id from filename
  sim_id <- as.numeric(gsub(".*sim_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Load dataset
  data <- readRDS(file_path)
  
  # Run VarSelCluster
  K <- data$params$n_clust
  vslcm <- VarSelCluster(x = data$X, gvals = K)
  
  # Calculate AMI and ARI
  ami_score <- AMI(data$labels, vslcm@partitions@zMAP)
  ari_score <- ARI(data$labels, vslcm@partitions@zMAP)
  
  # Save results
  result <- list(
    sim_id = sim_id,
    ami = ami_score,
    ari = ari_score,
    params = data$params,
    n_nonzero = length(vslcm@model@names.relevant)
  )
  
  output_file <- sprintf("varselcluster_res/varselcluster_%04d.rds", sim_id)
  saveRDS(result, output_file)
  
  return(sim_id)
}

# Process all datasets in parallel
n_cores <- 100
cat("Processing", length(dataset_files), "datasets using", n_cores, "cores...\n")

cl <- makeCluster(n_cores)
clusterExport(cl, c("process_dataset"))
clusterEvalQ(cl, {
  library(VarSelLCM)
  library(aricode)
})

results <- parLapply(cl, dataset_files, process_dataset)

stopCluster(cl)

cat("\nVarSelCluster analysis complete!\n")
cat("Results saved to: varselcluster_res/\n")