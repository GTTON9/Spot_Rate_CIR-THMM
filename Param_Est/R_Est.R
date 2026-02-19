estimate_R_sample_from_BW <- function(
    X,                       # d x (TT+1) or d x TT
    params_list,             # list(X1_BW_res$params, X2_BW_res$params, X3_BW_res$params)
    states_list = NULL,      # optional: list(X1_BW_res$states_estimate, ...)
    dt = 1/252,
    eps = 1e-8
) {
  stopifnot(is.matrix(X))
  d <- nrow(X)
  stopifnot(length(params_list) == d)
  
  TT <- ncol(X)
  has_increments <- (TT >= 2)
  
  # if X is d x (N+1): use increments; if it's d x N: still use increments
  dX   <- X[, 2:TT, drop = FALSE] - X[, 1:(TT-1), drop = FALSE]
  Xlag <- X[, 1:(TT-1), drop = FALSE]
  Tinc <- ncol(dX)
  
  # helper: pull (kappa, theta, sigma) for coord j at time k based on state
  get_pars_jk <- function(j, k) {
    p <- params_list[[j]]
    
    # default: use state 1 if no states_list supplied
    s <- 1L
    if (!is.null(states_list)) {
      svec <- states_list[[j]]
      # svec might be length Tinc or TT; align to increments
      s <- svec[min(k, length(svec))]
      # make sure it's 1..nstates
      s <- as.integer(s)
      if (s < 1) s <- 1
    }
    
    list(
      kappa = as.numeric(p$kappa[s]),
      theta = as.numeric(p$theta[s]),
      sigma = as.numeric(p$sigma[s])
    )
  }
  
  # standardized innovations E[k, j]
  E <- matrix(NA_real_, nrow = Tinc, ncol = d)
  
  for (k in 1:Tinc) {
    for (j in 1:d) {
      pars <- get_pars_jk(j, k)
      mu   <- pars$kappa * (pars$theta - Xlag[j, k]) * dt
      den  <- pars$sigma * sqrt(pmax(Xlag[j, k], eps)) * sqrt(dt)
      
      E[k, j] <- (dX[j, k] - mu) / den
    }
  }
  
  R_hat <- cor(E, use = "pairwise.complete.obs")
  R_hat <- (R_hat + t(R_hat)) / 2
  diag(R_hat) <- 1
  colnames(R_hat) <- rownames(R_hat) <- paste0("X", 1:d)
  
  list(R_hat = R_hat, innovations = E)
}