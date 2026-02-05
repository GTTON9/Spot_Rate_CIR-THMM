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



# 
# simulate_heston_dependent <- function(S0, v0, Reg_series, regime_params,
#                                    T, N, M = 1,
#                                    method = "E",
#                                    interp = TRUE, substeps = 100,
#                                    seed = 999, min_var = 1e-6) 