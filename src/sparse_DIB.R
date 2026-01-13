sparse_DIB <- function(X, ncl, s = -1, contkernel = "gaussian",
                       nomkernel = "aitchisonaitken",
                       ordkernel = "liracine",
                       cat_first = FALSE,
                       nstart = 100, maxiter = 100, maxsteps = 10,
                       randinit = NULL, verbose = FALSE,
                       sparsity, init_weights = "uniform",
                       n_cores = NULL){
  
  if (is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  if (init_weights == "warm"){
    km_test <- kmeans(x = X, centers = ncl, nstart = 1000)
  }
  
  X <- scale(X)
  contcols <- c(1:ncol(X))
  catcols <- c()
  tryCatch({
    suppressWarnings({
      s <- IBclust:::compute_bandwidth_cont(X, contkernel = contkernel)
      bws_vec <- rep(unique(s), ncol(X))
    })
  }, error = function(e) {
    bws_vec <<- rep(80, ncol(X))
  }, warning = function(w) {
    bws_vec <<- rep(80, ncol(X))
  })
  weights <- rep(1/sqrt(ncol(X)), ncol(X))
  converged <- FALSE
  step_count <- 0
  old_clust <- rep(-1, nrow(X))
  
  while (!converged && step_count < maxsteps){
    pxy_list <- coord_to_pxy_weights(X, s = bws_vec,
                                     weights = weights,
                                     contkernel = contkernel,
                                     use_parallel = FALSE,
                                     n_cores = n_cores,
                                     cluster = cl)
    py_x <- pxy_list$py_x
    px <- pxy_list$px
    pxy <- pxy_list$pxy
    hy <- pxy_list$hy
    
    if (step_count == 0 & init_weights == "warm"){
      best_clust <- list(Cluster = km_test$cluster)
    } else {
      best_clust <- IBclust:::DIBmix_iterate(X, ncl = ncl, randinit = best_clust$Cluster,
                                             tol = 0, py_x, hy, px, maxiter,
                                             bws_vec, contcols, catcols,
                                             runs = nstart, verbose = verbose)
    }
    cluster_labels <- best_clust$Cluster
    n <- length(best_clust$Cluster)
    K <- max(best_clust$Cluster)
    qt_x <- matrix(0, nrow = K, ncol = n)
    qt_x[cbind(best_clust$Cluster, 1:n)] <- 1
    qt_list <- IBclust:::qt_step(X, qt_x, ptol = 0, quiet = FALSE)
    qt <- qt_list$qt
    qt_x <- qt_list$qt_x
    
    s_val <- 1
    n_vars <- ncol(X)
    
    py_x_list <- lapply(1:n_vars, function(j) {
      compute_py_x_single(j, X, s_val, contkernel, nomkernel, ordkernel)
    })
    
    names(py_x_list) <- colnames(X)
    
    iyt_vec <- sapply(py_x_list, function(py_x) {
      qy_t <- IBclust:::qy_t_step_cpp(py_x, qt_x, qt, px)
      hy_t <- crossprod(qt, IBclust:::entropy(qy_t))
      hy - hy_t
    })
    
    old_weights <- weights
    weights <- dykstra_projection(iyt_vec, s = sparsity)
    weight_change <- sum(abs(weights - old_weights)) / sum(abs(old_weights))
    if(weight_change < 1e-5) { 
      converged <- TRUE
      break
    } else {
      old_weights <- weights
      old_clust <- best_clust$Cluster
    }
    step_count <- step_count + 1
  }
  return(list(best_clust = best_clust, weights = weights))
}