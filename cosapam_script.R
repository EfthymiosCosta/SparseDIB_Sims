library(aricode)
library(parallel)
library(rCOSA)
library(cluster)

# Create output directory
if (!dir.exists("cosapam_res")) dir.create("cosapam_res")

# Get list of all dataset files
dataset_files <- list.files("datasets", pattern = "^sim_\\d+\\.rds$", full.names = TRUE)
cat("Found", length(dataset_files), "datasets to process\n")

# Function to fix cosa2
fix_cosa2 <- function() {
  cosa2_fixed <- cosa2
  func_body <- deparse(body(cosa2))
  func_body <- gsub('class\\(X\\) == "data.frame"', 
                    'inherits(X, "data.frame")', 
                    func_body)
  body(cosa2_fixed) <- parse(text = paste(func_body, collapse = "\n"))[[1]]
  return(cosa2_fixed)
}

process_dataset <- function(file_path) {
  # Extract sim_id from filename
  sim_id <- as.numeric(gsub(".*sim_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Load dataset
  data <- readRDS(file_path)
  X <- as.data.frame(data$X)
  labels <- data$labels
  K <- data$params$n_clust
  
  # Fix cosa2 function
  cosa2_fixed <- fix_cosa2()
  
  # Grid search over lambda values
  lambda_seq <- seq(0.1, 2, by = 0.1)
  aris <- numeric(length(lambda_seq))
  amis <- numeric(length(lambda_seq))
  
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    
    # Run COSA
    cosa_rslts <- cosa2_fixed(X, lambda = lambda, targ = 'high/low', pwr = 2)
    Distmat <- as.matrix(cosa_rslts$D)
    
    # Run PAM clustering
    cosa_clust <- pam(Distmat, k = K, diss = TRUE, nstart = 100)
    
    # Calculate AMI and ARI
    aris[i] <- ARI(cosa_clust$clustering, labels)
    amis[i] <- AMI(cosa_clust$clustering, labels)
  }
  
  # Find best lambda
  best_idx <- which.max(aris)
  best_lambda <- lambda_seq[best_idx]
  best_ari <- aris[best_idx]
  best_ami <- amis[best_idx]
  
  # Save results
  result <- list(
    sim_id = sim_id,
    ari = best_ari,
    ami = best_ami,
    params = data$params
  )
  
  output_file <- sprintf("cosapam_res/cosapam_%04d.rds", sim_id)
  saveRDS(result, output_file)
  
  cat("Completed sim:", sim_id, '\n')
  
  return(sim_id)
}

n_cores <- 123
cat("Processing", length(dataset_files), "datasets using", n_cores, "cores...\n")

cl <- makeCluster(n_cores)
clusterExport(cl, c("process_dataset", "fix_cosa2"))
clusterEvalQ(cl, {
  library(aricode)
  library(rCOSA)
  library(cluster)
})

results <- parLapply(cl, dataset_files, process_dataset)
stopCluster(cl)

cat("\nCOSA+PAM analysis complete!\n")
cat("Results saved to: cosapam_res/\n")