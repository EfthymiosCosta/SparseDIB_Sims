library(RPEClust)

# Optimised version - initial one in RPEClust package has some bugs
fix_RPGMMClu_optimized <- function() {
  orig_func <- RPEClust::RPGMMClu
  env <- environment(orig_func)
  
  new_func <- function (x, true.cl = NULL, g, d = NULL, c = 10, B = 1000, B.star = 100, 
                        modelNames = NULL, diagonal = FALSE, ensmethod = "DWH", seed = 101, 
                        verb = FALSE) 
  {
    # Convert input to matrix if it's a data.frame
    x <- as.matrix(x)
    
    p <- ncol(x)
    if (is.null(d)) 
      d <- ceiling(c * log(g))
    else d = d
    n <- nrow(x)
    set.seed(seed)
    RPbase <- RPEClust:::generateRP(p = ncol(x), d = d, B = B)
    index = matrix(1:(d * B), d, B)
    Bic <- Ari <- bic.c <- bic.nc <- loglik.c <- det.s.ybar.y <- NULL
    cl.m <- matrix(NA, n, B)
    
    for (b in 1:B) {
      A <- as.matrix(RPbase[, index[, b]])
      A.bar <- qr.Q(qr(A), complete = TRUE)[, (d + 1):p]
      y <- x %*% A
      y.bar <- x %*% A.bar
      y.star <- cbind(y, y.bar)
      out <- mclust::Mclust(y, g, verbose = FALSE, modelNames = modelNames)
      loglik.c[b] <- out$loglik
      cl.m[, b] <- out$cl
      bic.c[b] <- out$bic
      X <- cbind(matrix(1, n), y)
      Y <- y.bar
      
      qr_X <- qr(X)
      B_mat <- qr.coef(qr_X, Y)
      
      Q <- qr.Q(qr_X)
      Q1 <- Q[, 1:qr_X$rank, drop = FALSE]
      residuals <- Y - Q1 %*% crossprod(Q1, Y)
      s.ybar.y <- crossprod(residuals) / (n - 1)
      
      if (diagonal) 
        s.ybar.y <- diag(diag(s.ybar.y))
      
      m <- ncol(y.bar)
      if (diagonal) 
        det.o <- prod(diag(s.ybar.y))
      else det.o <- det(s.ybar.y)
      
      if ((det.o) < 9.88131291682493e-324) 
        det.s.ybar.y[b] <- 9.88131291682493e-324
      else if (det.o == Inf) 
        det.s.ybar.y[b] <- 1e+308
      else det.s.ybar.y[b] <- det.o
      
      tryCatch({
        chol_s <- chol(s.ybar.y)
        z <- backsolve(chol_s, t(residuals), transpose = TRUE)
        quad_form <- sum(z^2)
      }, error = function(e) {
        # Fallback if Cholesky fails
        quad_form <<- sum(diag(residuals %*% solve(s.ybar.y, tol = 1e-30) %*% t(residuals)))
      })
      
      loglik.nc <- (-(m * n)/2) * log((2 * pi)) - (n/2) * log(det.s.ybar.y[b]) - 0.5 * quad_form
      
      if (diagonal) 
        k <- m * (d + 1) + m
      else k <- m * (d + 1) + (m * (m + 1))/2
      
      bic.nc[b] <- 2 * loglik.nc - k * log(n)
      Bic[b] <- bic.c[b] + bic.nc[b]
      
      if (!is.null(true.cl)) 
        Ari[b] = mclust::adjustedRandIndex(out$classification, true.cl)
      
      if (verb) 
        print(b)
    }
    
    cl.ens.1 <- data.frame(cl.m[, order(Bic, decreasing = TRUE)[1:B.star]])
    cl.ens.1.2 <- lapply(cl.ens.1, function(x) clue::as.cl_membership(x))
    cl.consensus <- apply(clue::cl_consensus(cl.ens.1.2, method = ensmethod)$.Data, 
                          1, which.max)
    if (!is.null(true.cl)) 
      ari <- mclust::adjustedRandIndex(cl.consensus, true.cl)
    else
      ari <- NULL
    
    if (!is.null(ari))
      names(ari) <- paste0("B.star=", B.star)
    
    ensemble <- list(ari = ari, label.vec = cl.consensus)
    individual <- list(label.vec = cl.m, ari = Ari, bic = Bic, 
                       bic.GMM = bic.c, bic.reg = bic.nc)
    return(list(ensemble = ensemble, individual = individual))
  }
  
  environment(new_func) <- env
  return(new_func)
}

RPGMMClu <- fix_RPGMMClu_optimized()

# Create output directory
if (!dir.exists("data_res")) dir.create("data_res")

bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.matrix(bladder_cancer$X)
labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"

K <- length(unique(labels_aggr))
cat("Running RPEClust with K =", K, "clusters...\n")

# Wrapper function with retry logic
run_rpeclust_robust <- function(x, true.cl, g, B, max_attempts = 10) {
  for (attempt in 1:max_attempts) {
    seed_val <- attempt  
    cat("\nAttempt", attempt, "with seed", seed_val, "...\n")
    
    result <- tryCatch({
      RPGMMClu(x = x, true.cl = true.cl, g = g, B = B, seed = seed_val, verb = TRUE)
    }, error = function(e) {
      cat("  Error encountered:", e$message, "\n")
      return(NULL)
    }, warning = function(w) {
      cat("  Warning:", w$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      cat("  Success!\n")
      return(result)
    }
  }
  
  stop("Failed to complete RPGMMClu after ", max_attempts, " attempts")
}

cat("\nStarting clustering analysis...\n")
res <- run_rpeclust_robust(x = X, true.cl = labels_aggr, g = K, B = 500)  

output_file <- "data_res/rpeclust_bladder.RDS"
saveRDS(res$ensemble$label.vec, output_file)

cat("\nRPEClust analysis complete!\n")
cat("Clustering saved to:", output_file, "\n")