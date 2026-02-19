simulate_Regime <- function(Gamma, N, n_dim = 1, init_state = 0, do_plot = TRUE) {
  
  if (length(init_state) == 1) init_state <- rep(init_state, n_dim)
  stopifnot(length(init_state) == n_dim)
  
  Reg_mat <- matrix(0, nrow = n_dim, ncol = N)
  Reg_mat[, 1] <- init_state
  
  for (j in 1:n_dim) {
    for (t in 2:N) {
      Reg_mat[j, t] <- sample(
        c(0, 1),
        1,
        prob = Gamma[Reg_mat[j, t-1] + 1, ]
      )
    }
  }
  
  
  return(Reg_mat)
}





simulate_Regime_multi <- function(Gamma, N, init_state = 0) {
  # Gamma: list of length n_dim, each is KxK transition matrix
  stopifnot(is.list(Gamma), length(Gamma) >= 1)
  n_dim <- length(Gamma)
  
  K <- nrow(Gamma[[1]])
  stopifnot(ncol(Gamma[[1]]) == K)
  
  # sanity check all Gamma matrices
  for (j in 1:n_dim) {
    Gj <- as.matrix(Gamma[[j]])
    stopifnot(nrow(Gj) == K, ncol(Gj) == K)
    if (max(abs(rowSums(Gj) - 1)) > 1e-8) stop("Each Gamma[[j]] must have rows summing to 1.")
    if (min(Gj) < -1e-12) stop("Each Gamma[[j]] must be nonnegative (within tolerance).")
    Gamma[[j]] <- Gj
  }
  
  # init_state: scalar or length n_dim, values in 0..K-1
  if (length(init_state) == 1) init_state <- rep(init_state, n_dim)
  stopifnot(length(init_state) == n_dim)
  if (min(init_state) < 0 || max(init_state) > (K-1)) stop("init_state must be in 0..K-1.")
  
  Reg_mat <- matrix(0, nrow = n_dim, ncol = N)
  Reg_mat[, 1] <- init_state
  
  states <- 0:(K-1)
  
  for (j in 1:n_dim) {
    Gj <- Gamma[[j]]
    for (t in 2:N) {
      prev <- Reg_mat[j, t-1]
      Reg_mat[j, t] <- sample(states, 1, prob = Gj[prev + 1, ])
    }
  }
  
  Reg_mat
}





simulate_heston_simple <- function(S0, v0, Reg_series, regime_params,
                            T, N, M = 1,
                            method = "E",
                            interp = TRUE, substeps = 100,
                            seed = 999, min_var = 1e-6) {
  
  if (!is.null(seed)) set.seed(seed)
  method <- toupper(method)
  
  # ----------------------------
  # check regime_params format
  # ----------------------------
  stopifnot(is.list(regime_params), length(regime_params) >= 2)
  
  need_names <- c("mu","kappa","theta","sigma","rho")
  for (rp in regime_params) {
    if (!all(need_names %in% names(rp)))
      stop("Each regime_params element must contain: mu,kappa,theta,sigma,rho")
  }
  
  # Feller warning
  for (r in seq_along(regime_params)) {
    p <- regime_params[[r]]
    if (2*p$kappa*p$theta < p$sigma^2) {
      message(paste("Warning: Feller condition violated in regime", r))
    }
  }
  
  # ----------------------------
  # accept vector or matrix regime
  # ----------------------------
  if (is.vector(Reg_series)) {
    Reg_series <- matrix(Reg_series, nrow = 1)
  } else {
    Reg_series <- as.matrix(Reg_series)
  }
  
  n_dim <- nrow(Reg_series)
  stopifnot(ncol(Reg_series) == N)
  
  # ----------------------------
  # interpolation
  # ----------------------------
  if (interp) {
    N_use <- N * substeps
    Reg_use <- Reg_series[, rep(1:N, each=substeps), drop=FALSE]
  } else {
    N_use <- N
    Reg_use <- Reg_series
  }
  
  dt <- T / N_use
  sqrt_dt <- sqrt(dt)
  
  # ----------------------------
  # output containers
  # ----------------------------
  if (M == 1) {
    S_out <- matrix(0, n_dim, N_use+1)
    V_out <- matrix(0, n_dim, N_use+1)
    S_out[,1] <- S0
    V_out[,1] <- v0
  } else {
    S_out <- array(0, c(n_dim, M, N_use+1))
    V_out <- array(0, c(n_dim, M, N_use+1))
    S_out[,,1] <- S0
    V_out[,,1] <- v0
  }
  
  # ============================
  # simulate per dimension
  # ============================
  for (j in 1:n_dim) {
    
    if (M == 1) {
      S_path <- numeric(N_use+1); S_path[1] <- S0
      V_path <- numeric(N_use+1); V_path[1] <- v0
    } else {
      S_path <- matrix(0, M, N_use+1); S_path[,1] <- S0
      V_path <- matrix(0, M, N_use+1); V_path[,1] <- v0
    }
    
    for (i in 1:N_use) {
      
      reg_id <- Reg_use[j,i] + 1
      p <- regime_params[[reg_id]]
      
      mu    <- p$mu
      kappa <- p$kappa
      theta <- p$theta
      sigma <- p$sigma
      rho   <- p$rho
      
      rho_p <- sqrt(max(0, 1-rho^2))
      
      Z1 <- rnorm(M)
      Z2 <- rnorm(M)
      
      dW1 <- Z1 * sqrt_dt
      dW2 <- (rho*Z1 + rho_p*Z2) * sqrt_dt
      
      V_prev <- if (M==1) V_path[i] else V_path[,i]
      V_pos  <- pmax(V_prev, 0)
      sqrtV  <- sqrt(V_pos)
      
      # ---- variance update ----
      if (method == "E_C") {
        
        exp_k_dt <- exp(-kappa*dt)
        C_val <- (2*kappa)/((1-exp_k_dt)*(sigma^2))
        df_val <- 4*kappa*theta/sigma^2
        lambda <- 2*C_val*V_pos*exp_k_dt
        
        Y <- rchisq(M, df=df_val, ncp=lambda)
        V_next <- Y/(2*C_val)
        
      } else if (method == "E") {
        
        V_next <- V_prev +
          kappa*(theta-V_prev)*dt +
          sigma*sqrtV*dW2
        
      } else if (method == "M") {
        
        V_e <- V_prev +
          kappa*(theta-V_prev)*dt +
          sigma*sqrtV*dW2
        
        V_next <- V_e + 0.25*sigma^2*(dW2^2-dt)
        
      } else stop("method must be E, M, or E_C")
      
      V_next <- pmax(V_next, min_var)
      
      # ---- price update ----
      S_prev <- if (M==1) S_path[i] else S_path[,i]
      
      S_next <- S_prev * exp(
        (mu - 0.5*V_pos)*dt +
          sqrtV*dW1
      )
      
      if (M==1) {
        V_path[i+1] <- V_next
        S_path[i+1] <- S_next
      } else {
        V_path[,i+1] <- V_next
        S_path[,i+1] <- S_next
      }
    }
    
    if (M==1) {
      S_out[j,] <- S_path
      V_out[j,] <- V_path
    } else {
      S_out[j,,] <- S_path
      V_out[j,,] <- V_path
    }
  }
  
  return(list(
    S_paths = S_out,
    V_paths = V_out,
    dt = dt,
    method = method
  ))
}












simulate_CIR_dependent <- function(x0, Reg_series, param_mat,
                                      T, N, M = 1,
                                      interp = TRUE, substeps = 100,
                                      seed = 999, min_state = 1e-8,
                                      R = NULL,
                                      alpha0 = 0.001,
                                      alpha = NULL) {
  # x0: length d initial states
  # Reg_series: d x N (values 0..K-1). Here assume K=2, values 0/1.
  # param_mat: (d*K) x 3 matrix with columns kappa, theta, sigma
  #            rows: regime1 for all d, then regime2 for all d (block by regime)
  # R: d x d correlation matrix for Brownian increments
  # alpha0, alpha: short-rate mapping r_t = alpha0 + alpha' X_t
  # alpha: length d, default 0 (or user supplied)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Reg_series -> matrix d x N
  if (is.vector(Reg_series)) Reg_series <- matrix(Reg_series, nrow = 1)
  Reg_series <- as.matrix(Reg_series)
  d <- nrow(Reg_series)
  stopifnot(ncol(Reg_series) == N)
  if (length(x0) != d) stop("x0 must have length d = nrow(Reg_series).")

  
  param_mat <- as.matrix(param_mat)
 
   colnames(param_mat) <- colnames(param_mat) %||% c("kappa","theta","sigma")  # keep if set
   
  
  # correlation matrix
  if (is.null(R)) {
    R <- diag(d)
  } else {
    R <- as.matrix(R)
    stopifnot(nrow(R) == d, ncol(R) == d)
    if (max(abs(diag(R) - 1)) > 1e-8) stop("R must have diag=1.")
    chol(R)  # PD check
  }
  L <- chol(R)  
  
  # interpolation
  if (interp) {
    N_use <- N * substeps
    Reg_use <- Reg_series[, rep(1:N, each = substeps), drop = FALSE]
  } else {
    N_use <- N
    Reg_use <- Reg_series
  }
  
  dt <- T / N_use
  sqrt_dt <- sqrt(dt)
  
  # output X
  if (M == 1) {
    X_out <- matrix(0, d, N_use + 1)
    X_out[, 1] <- x0
  } else {
    X_out <- array(0, c(d, M, N_use + 1))
    for (j in 1:d) X_out[j, , 1] <- x0[j]
  }
  
  # output r_t
  if (M == 1) {
    r_out <- numeric(N_use + 1)
    r_out[1] <- alpha0 + sum(alpha * x0)
  } else {
    r_out <- matrix(0, M, N_use + 1)
    r_out[, 1] <- alpha0 + as.numeric(alpha %*% matrix(x0, ncol = 1))
  }
  
  # helper to get row in param_mat for (regime, series)
  # regime index in code: 1..K
  idx_row <- function(reg_id, j) {
    (reg_id - 1) * d + j
  }
  
  # simulate Euler with correlated Brownian drivers
  for (i in 1:N_use) {
    Z_unc  <- matrix(rnorm(M * d), nrow = M, ncol = d)
    Z_corr <- Z_unc %*% L
    dW_mat <- Z_corr * sqrt_dt
    
    for (j in 1:d) {
      reg_id <- Reg_use[j, i] + 1  # 1..K
      row_id <- idx_row(reg_id, j)
      
      kappa <- param_mat[row_id, 1]
      theta <- param_mat[row_id, 2]
      sigma <- param_mat[row_id, 3]
      
      X_prev <- if (M == 1) X_out[j, i] else X_out[j, , i]
      X_pos  <- pmax(X_prev, 0)
      sqrtX  <- sqrt(X_pos)
      
      dWj <- if (M == 1) dW_mat[1, j] else dW_mat[, j]
      
      X_next <- X_prev +
        kappa * (theta - X_prev) * dt +
        sigma * sqrtX * dWj
      
      X_next <- pmax(X_next, min_state)
      
      if (M == 1) {
        X_out[j, i + 1] <- X_next
      } else {
        X_out[j, , i + 1] <- X_next
      }
    }
    
    # update r_t at i+1 using X_{i+1}
    if (M == 1) {
      r_out[i + 1] <- alpha0 + sum(alpha * X_out[, i + 1])
    } else {
      # X_out[, , i+1] is d x M; want M-vector of alpha0 + alpha'X
      X_now <- X_out[, , i + 1, drop = FALSE]  # d x M x 1
      X_now <- X_now[, , 1]                    # d x M
      r_out[, i + 1] <- alpha0 + as.numeric(t(alpha) %*% X_now)
    }
  }
  
  list(
    X_paths = X_out,
    r_paths = r_out,
    dt = dt,
    corr = R,
    interp = interp,
    substeps = substeps
  )
}


`%||%` <- function(a, b) if (!is.null(a)) a else b



