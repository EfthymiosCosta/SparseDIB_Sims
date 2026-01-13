library(mclust)
library(parallel)
library(sparcl)
if (!dir.exists("skmeans_res")) dir.create("skmeans_res")
dataset_files <- list.files("datasets", pattern = "^sim_\\d+\\.rds$", full.names = TRUE)
cat("Found", length(dataset_files), "datasets to process\n")
process_dataset <- function(file_path) {
  sim_id <- as.numeric(gsub(".*sim_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  data <- readRDS(file_path)
  X <- data$X
  labels <- data$labels
  ncl <- data$params$n_clust
  
  wbound_seq <- seq(1.2, 10, by = 0.1)
  aris <- numeric(length(wbound_seq))
  amis <- numeric(length(wbound_seq))
  n_nonzero <- numeric(length(wbound_seq))
  
  for (i in seq_along(wbound_seq)) {
    wbound <- wbound_seq[i]
    
    result <- tryCatch({
      SKM <- KMeansSparseCluster(X, K = ncl, wbounds = wbound, maxiter = 100,
                                 nstart = 100, silent = TRUE)
      
      list(
        ari = adjustedRandIndex(SKM[[1]]$Cs, labels),
        ami = aricode::AMI(SKM[[1]]$Cs, labels),
        n_nonzero = sum(SKM[[1]]$ws > 0),
        success = TRUE
      )
    }, error = function(e) {
      cat("Warning: wbound =", wbound, "failed for sim", sim_id, "\n")
      list(
        ari = NA,
        ami = NA,
        n_nonzero = NA,
        success = FALSE
      )
    })
    
    if (result$success) {
      aris[i] <- result$ari
      amis[i] <- result$ami
      n_nonzero[i] <- result$n_nonzero
    } else {
      aris[i] <- NA
      amis[i] <- NA
      n_nonzero[i] <- NA
    }
  }
  
  result <- list(
    sim_id = sim_id,
    ari = aris,
    ami = amis,
    wbound = wbound_seq,
    n_nonzero = n_nonzero,
    params = data$params
  )
  
  output_file <- sprintf("skmeans_res/skmeans_%04d.rds", sim_id)
  saveRDS(result, output_file)
  cat("Completed sim", sim_id, "\n")
  
  return(sim_id)
}
n_cores <- 123
cat("Processing", length(dataset_files), "datasets using", n_cores, "cores...\n")
cl <- makeCluster(n_cores)
clusterExport(cl, c("process_dataset"))
clusterEvalQ(cl, {
  library(mclust)
  library(sparcl)
  library(aricode)
})
results <- parLapply(cl, dataset_files, process_dataset)
stopCluster(cl)
cat("\nSparse K-Means analysis complete!\n")
cat("Results saved to: skmeans_res/\n")