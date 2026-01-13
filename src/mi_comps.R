library(parallel)

# Compute marginal KDEs
compute_marginals <- function(X, parallel = TRUE, n_cores = NULL) {
  p <- ncol(X)
  n <- nrow(X)
  
  compute_single_marginal <- function(i) {
    data_i <- data.frame(X = X[, i])
    bw <- np::npudensbw(~X, data = data_i)
    f <- np::npksum(txdat = X[, i],
                exdat = X[, i],
                bws = bw$bw,
                bandwidth.divide = TRUE)$ksum / n
    list(f = f, bw = bw$bw)
  }
  
  if (parallel) {
    if (is.null(n_cores)) n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("X", "n"), envir = environment())
    clusterEvalQ(cl, library(np))
    marginals <- parLapply(cl, 1:p, compute_single_marginal)
    stopCluster(cl)
  } else {
    marginals <- lapply(1:p, compute_single_marginal)
  }
  return(marginals)
}

# MI function using precomputed marginals
compute_mi_fast <- function(x1, x2, f_x1, f_x2) {
  data <- data.frame(X1 = x1, X2 = x2)
  n <- nrow(data)
  
  # Estimate joint density f(X1, X2)
  bw_joint <- np::npudensbw(~X1 + X2, data = data)
  f_joint <- np::npksum(txdat = data, 
                    exdat = data,
                    bws = bw_joint$bw,
                    bandwidth.divide = TRUE)$ksum / n
  
  # Compute MI
  epsilon <- 1e-10
  mi <- mean(log((f_joint + epsilon) / ((f_x1 + epsilon) * (f_x2 + epsilon))))
  
  return(mi)
}

# Main function to compute pairwise MIs
compute_pairwise_mi <- function(X, parallel = TRUE, n_cores = NULL) {
  p <- ncol(X)
  var_names <- colnames(X)
  if (is.null(var_names)) {
    var_names <- paste0("V", 1:p)
  }
  
  if (parallel){
    marginals <- compute_marginals(X, parallel = TRUE, n_cores = n_cores)
  } else {
    marginals <- compute_marginals(X, parallel = FALSE)
  }
  
  pairs <- combn(p, 2)
  n_pairs <- ncol(pairs)
  
  compute_pair_mi <- function(pair_idx) {
    i <- pairs[1, pair_idx]
    j <- pairs[2, pair_idx]
    
    mi_val <- tryCatch({
      compute_mi_fast(X[, i], X[, j], 
                      marginals[[i]]$f, 
                      marginals[[j]]$f)
    }, error = function(e) {
      warning(paste("Error computing MI for pair (", i, ",", j, "):", e$message))
      NA
    })
    
    return(list(i = i, j = j, mi = mi_val))
  }
  
  # Compute in parallel or sequential
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- detectCores() - 1
    }
    cat("Computing", n_pairs, "pairwise MIs using", n_cores, "cores...\n")
    
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("compute_mi_fast", "X", "marginals", "pairs"), 
                  envir = environment())
    clusterEvalQ(cl, library(np))
    
    results <- parLapply(cl, 1:n_pairs, compute_pair_mi)
    stopCluster(cl)
  } else {
    cat("Computing", n_pairs, "pairwise MIs sequentially...\n")
    results <- lapply(1:n_pairs, compute_pair_mi)
  }
  
  # Build MI matrix
  mi_matrix <- matrix(0, nrow = p, ncol = p)
  rownames(mi_matrix) <- var_names
  colnames(mi_matrix) <- var_names
  
  for (result in results) {
    if (!is.na(result$mi)) {
      mi_matrix[result$i, result$j] <- result$mi
      mi_matrix[result$j, result$i] <- result$mi
    }
  }
  
  # Set negative values to zero (caused by numeric overflow)
  mi_matrix[mi_matrix < 0] <- 0
  
  return(mi_matrix)
}


# Function to compute py_x for a single variable
compute_py_x_single <- function(j, X, s_val, contkernel, nomkernel, ordkernel) {
  X_j <- as.data.frame(X[, j, drop = FALSE])
  # Compute pxy
  pxy_list <- IBclust:::coord_to_pxy_R(X_j,
                                       s = s_val,
                                       lambda = -1,
                                       cat_cols = integer(0),
                                       cont_cols = 1,
                                       contkernel = contkernel,
                                       nomkernel = nomkernel,
                                       ordkernel = ordkernel)
  
  return(pxy_list$py_x)
}

# Function to compute iyt for a single py_x
compute_iyt_single <- function(py_x, qt_x, qt, px, hy) {
  qy_t <- IBclust:::qy_t_step_cpp(py_x, qt_x, qt, px)
  ht <- IBclust:::entropy(qt)
  hy_t <- crossprod(qt, IBclust:::entropy(qy_t))
  iyt <- hy - hy_t
  
  return(iyt)
}
