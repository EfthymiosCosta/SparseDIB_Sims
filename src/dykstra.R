# Project onto nonnegative L1-ball using Duchi's algorithm
projection_nonnegative_l1_ball <- function(x, s) {
  x <- pmax(x, 0)
  if (sum(x) <= s) {
    return(x)
  }
  # Project onto simplex using sorting
  u <- sort(x, decreasing = TRUE)
  cssv <- cumsum(u) - s
  ind <- 1:length(x)
  cond <- u - cssv / ind > 0
  rho <- max(ind[cond])
  theta <- cssv[rho] / rho
  return(pmax(x - theta, 0))
}

# Project onto nonnegative L2-ball
projection_nonnegative_l2_ball <- function(x) {
  x <- pmax(x, 0)
  norm_val <- norm(x, type = "2")
  if (norm_val <= 1) {
    return(x)
  }
  return(x / norm_val)
}

dykstra_projection <- function(y, s, max_iter = 1000, tol = 1e-8) {
  n <- length(y)
  x <- y
  # L1 and L2 projection residuals
  p <- numeric(n)
  q <- numeric(n)
  
  for (i in 1:max_iter) {
    x_prev <- x
    
    # L1 projection
    y1 <- projection_nonnegative_l1_ball(x + p, s)
    p <- x + p - y1
    
    # L2 projection  
    x <- projection_nonnegative_l2_ball(y1 + q)
    q <- y1 + q - x
    
    # Check convergence
    if (norm(x - x_prev, type = "2") < tol) {
      break
    }
  }
  
  return(x)
}
