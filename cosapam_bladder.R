library(aricode)
library(parallel)
library(rCOSA)
library(cluster)

# Create output directory
if (!dir.exists("data_res")) dir.create("data_res")

# Bladder cancer dataset
bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.data.frame(bladder_cancer$X)
dim(X)
labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"

# Determine number of clusters from labels
K <- length(unique(labels_aggr))

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

# Lambda sequence from 0.1 to 5.0
lambda_seq <- seq(0.1, 5, by = 0.1)

# Function to process a single lambda value
process_lambda <- function(lambda) {
  # Fix cosa2 function
  cosa2_fixed <- fix_cosa2()
  
  # Run COSA
  cosa_rslts <- cosa2_fixed(X, lambda = lambda, targ = 'high/low', pwr = 2)
  Distmat <- as.matrix(cosa_rslts$D)
  
  # Run PAM clustering
  cosa_clust <- pam(Distmat, k = K, diss = TRUE, nstart = 100)
  
  # Save only the clustering vector
  output_file <- sprintf("data_res/cosapam_bladder_%.1f.RDS", lambda)
  saveRDS(cosa_clust$clustering, output_file)
  
  cat("Completed lambda:", lambda, '\n')
  
  return(lambda)
}

# Parallel processing
n_cores <- 100
cl <- makeCluster(n_cores)
clusterExport(cl, c("X", "K", "fix_cosa2"))
clusterEvalQ(cl, {
  library(rCOSA)
  library(cluster)
})

results <- parLapply(cl, lambda_seq, process_lambda)
stopCluster(cl)

cat("All lambda values completed!\n")