source('ALYZ.R')
library(MASS)

if (!dir.exists("datasets")) dir.create("datasets")

set.seed(1234) 
n <- 200
n_cols <- c(0.5, 1, 2, 5) * n
ratio_relevant <- c(0.05, 0.10, 0.20, 0.50)
n_clust <- c(3, 5, 8)
sph <- c(TRUE, FALSE)
balanced <- c(TRUE, FALSE)
n_replicates <- 50

# Parameter combinations
param_grid <- expand.grid(
  n = n,
  p = n_cols,
  ratio_relevant = ratio_relevant,
  n_clust = n_clust,
  spherical = sph,
  balanced = balanced,
  replicate = 1:n_replicates
)

# Simulation ID and seed
param_grid$sim_id <- 1:nrow(param_grid)
param_grid$seed <- sample.int(.Machine$integer.max, nrow(param_grid))

param_grid <- param_grid[, c("sim_id", "n", "p", "ratio_relevant", 
                             "n_clust", "spherical", "balanced", 
                             "replicate", "seed")]

write.csv(param_grid, "datasets/simulation_index.csv", row.names = FALSE)

generate_cluster_means <- function(K, q, spherical) {
  # Same separation for both spherical and non-spherical
  if (K == 3) {
    multipliers <- c(-1, 0, 1)
  } else if (K == 5) {
    multipliers <- c(-2, -1, 0, 1, 2)
  } else if (K == 8) {
    multipliers <- seq(-3.5, 3.5, by = 1)
  }
  
  mu_base <- 1.2
  
  means <- matrix(0, nrow = K, ncol = q)
  for (k in 1:K) {
    means[k, ] <- multipliers[k] * mu_base
  }
  
  return(means)
}

generate_dataset <- function(n, p, ratio_relevant, n_clust, spherical, balanced, seed) {
  set.seed(seed)
  
  q <- round(p * ratio_relevant)  
  
  if (balanced) {
    cluster_sizes <- rep(n %/% n_clust, n_clust)
    remainder <- n %% n_clust
    if (remainder > 0) {
      cluster_sizes[1:remainder] <- cluster_sizes[1:remainder] + 1
    }
  } else {
    min_size <- 10
    remaining <- n - (n_clust * min_size)
    
    if (remaining < 0) {
      stop("Cannot create ", n_clust, " clusters with minimum size ", 
           min_size, " and n = ", n)
    }
    cluster_sizes <- rep(min_size, n_clust)
    if (remaining > 0) {
      props <- rgamma(n_clust, shape = 0.5, rate = 1)
      props <- props / sum(props)
      additional <- round(remaining * props)
      cluster_sizes <- cluster_sizes + additional
      diff <- n - sum(cluster_sizes)
      if (diff != 0) {
        cluster_sizes[which.max(cluster_sizes)] <- cluster_sizes[which.max(cluster_sizes)] + diff
      }
    }
  }
  
  cluster_means <- generate_cluster_means(n_clust, q, spherical)
  X <- matrix(0, nrow = n, ncol = p)
  labels <- rep(1:n_clust, times = cluster_sizes)
  
  if (spherical) {
    for (k in 1:n_clust) {
      idx <- which(labels == k)
      n_k <- length(idx)
      for (j in 1:q) {
        X[idx, j] <- rnorm(n_k, mean = cluster_means[k, j], sd = 1)
      }
      if (q < p) {
        X[idx, (q+1):p] <- matrix(rnorm(n_k * (p - q), mean = 0, sd = 1), 
                                  nrow = n_k, ncol = p - q)
      }
    }
  } else {
    covariance_matrices <- vector("list", n_clust)
    for (k in 1:n_clust) {
      covariance_matrices[[k]] <- generate_covmat_ALYZ(q, CN = 200)
    }
    
    for (k in 1:n_clust) {
      idx <- which(labels == k)
      n_k <- length(idx)
      mean_vec <- cluster_means[k, ]
      X[idx, 1:q] <- mvrnorm(n = n_k, 
                             mu = mean_vec, 
                             Sigma = covariance_matrices[[k]])
      if (q < p) {
        X[idx, (q+1):p] <- matrix(rnorm(n_k * (p - q), mean = 0, sd = 1), 
                                  nrow = n_k, ncol = p - q)
      }
    }
  }
  
  list(
    X = X,
    labels = labels,
    params = list(
      n = n,
      p = p,
      q = q,
      ratio_relevant = ratio_relevant,
      n_clust = n_clust,
      spherical = spherical,
      balanced = balanced,
      cluster_sizes = cluster_sizes
    )
  )
}

# Generate all datasets in parallel
library(parallel)
n_cores <- 100
cat("Generating", nrow(param_grid), "datasets using", n_cores, "cores...\n")

generate_and_save <- function(i, param_grid) {
  row <- param_grid[i, ]
  
  dataset <- generate_dataset(
    n = row$n,
    p = row$p,
    ratio_relevant = row$ratio_relevant,
    n_clust = row$n_clust,
    spherical = row$spherical,
    balanced = row$balanced,
    seed = row$seed
  )
  
  filename <- sprintf("datasets/sim_%04d.rds", row$sim_id)
  saveRDS(dataset, filename)
  
  return(row$sim_id)
}

cl <- makeCluster(n_cores)
clusterExport(cl, c("generate_dataset", "generate_cluster_means", 
                    "param_grid", "generate_covmat_ALYZ"))
clusterEvalQ(cl, library(MASS))

results <- parLapply(cl,
                     1:nrow(param_grid),
                     generate_and_save, 
                     param_grid = param_grid)

stopCluster(cl)

cat("\nDataset generation complete!\n")