# Computes perturbed similarity matrix with weights in exponents
coord_to_pxy_weights <- function(X, s, weights, contkernel = "gaussian",
                                 use_parallel = TRUE, n_cores = NULL,
                                 cluster = NULL){ 
  n_features <- ncol(X)
  n_obs <- nrow(X)
  if (is.null(n_cores)) {
    n_cores <- min(parallel::detectCores() - 1, n_features)
  }
  use_parallel <- use_parallel && n_cores > 1
  
  compute_feature_pyx <- function(j) {
    bw_j <- s[j]
    X_j <- X[, j]
    
    kw <- np::npksum(bws = bw_j,
                     txdat = X_j,
                     exdat = X_j,
                     ckertype = contkernel,
                     return.kernel.weights = TRUE)$kw
    
    py_x_j <- t(kw)
    py_x_j <- sweep(py_x_j, 2, colSums(py_x_j), '/')
    
    return(log(py_x_j + 1e-10)) 
  }
  
  if (use_parallel) {
    # Use provided cluster or create new one
    if (is.null(cluster)) {
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterEvalQ(cl, library(np))
      created_cluster <- TRUE
    } else {
      cl <- cluster
      # Load np on the existing cluster
      parallel::clusterEvalQ(cl, {
        if (!"np" %in% .packages()) library(np)
      })
      created_cluster <- FALSE
    }
    
    parallel::clusterExport(cl, c("X", "s", "contkernel"), envir = environment())
    
    log_py_x_list <- parallel::parLapply(cl, 1:n_features, compute_feature_pyx)
    
    # Only stop if we created it
    if (created_cluster) {
      parallel::stopCluster(cl)
    }
  } else {
    log_py_x_list <- lapply(1:n_features, compute_feature_pyx)
  }
  
  log_py_x_combined <- log_py_x_list[[1]] * weights[1]
  
  for (j in 2:n_features) {
    log_py_x_combined <- log_py_x_combined + weights[j] * log_py_x_list[[j]]
    if (j %% 100 == 0) gc() 
  }
  
  log_py_x_normalized <- sweep(log_py_x_combined, 2, 
                               apply(log_py_x_combined, 2, max), "-")
  
  py_x <- exp(log_py_x_normalized)
  zero_rows <- which(rowSums(py_x) == 0 | !is.finite(rowSums(py_x)))
  if (length(zero_rows) > 0) {
    py_x <- py_x[-zero_rows, , drop = FALSE]
  }
  py_x <- sweep(py_x, 2, colSums(py_x), '/')
  px <- matrix(1/n_obs, nrow = nrow(py_x), ncol = n_obs)
  pxy <- t(py_x * px)
  
  hx <- IBclust:::entropySingle(rowSums(pxy))
  hy <- IBclust:::entropySingle(colSums(pxy))
  hy_x <- rowSums(pxy) %*% IBclust:::entropy(py_x)
  ixy <- hy - hy_x
  px <- matrix(1/n_obs, nrow = n_obs, ncol = 1)
  
  return(list('py_x' = py_x, 'px' = px, 'pxy' = pxy, 'hy' = hy))
}