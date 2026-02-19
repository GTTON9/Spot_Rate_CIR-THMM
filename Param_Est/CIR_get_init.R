initial_heuristic_multi <- function(X, states = 2, dt = 1/252,
                                    order_by = c("theta", "sigma"),
                                    calm_first = TRUE) {
  stopifnot(is.matrix(X))
  order_by <- match.arg(order_by)
  
  d  <- nrow(X)
  TT <- ncol(X)
  
  cluster_mat <- matrix(NA_integer_, nrow = d, ncol = TT)
  Gamma_list  <- vector("list", d)
  
  kappa_mat <- matrix(NA_real_, nrow = states, ncol = d)
  theta_mat <- matrix(NA_real_, nrow = states, ncol = d)
  sigma_mat <- matrix(NA_real_, nrow = states, ncol = d)
  
  for (j in 1:d) {
    x <- as.numeric(X[j, ])
    
    ok <- !is.na(x)
    km <- stats::kmeans(x[ok], centers = states, iter.max = 100, nstart = 50)
    
    cluster <- rep(NA_integer_, TT)
    cluster[ok] <- km$cluster
    
    # fill NA in cluster by carry-forward
    if (anyNA(cluster)) {
      first_ok <- which(!is.na(cluster))[1]
      cluster[1:(first_ok - 1)] <- cluster[first_ok]
      if (first_ok < TT) {
        for (t in (first_ok + 1):TT) {
          if (is.na(cluster[t])) cluster[t] <- cluster[t - 1]
        }
      }
    }
    
    # --------- estimate theta/kappa/sigma (OLS) ----------
    kappa_j <- rep(NA_real_, states)
    theta_j <- rep(NA_real_, states)
    sigma_j <- rep(NA_real_, states)
    
    dx   <- diff(x)
    xlag <- x[-length(x)]
    clag <- cluster[-length(cluster)]
    
    for (s in 1:states) {
      theta_s <- mean(x[cluster == s], na.rm = TRUE)
      theta_j[s] <- max(theta_s, 1e-6)
      
      idx_s <- which(clag == s)
      if (length(idx_s) < 3) {
        kappa_j[s] <- 2.0
        sigma_j[s] <- 0.1
        next
      }
      
      dx_s   <- dx[idx_s]
      xlag_s <- xlag[idx_s]
      
      fit <- lm(dx_s ~ xlag_s)
      b   <- coef(fit)[["xlag_s"]]
      kappa_j[s] <- max(-b / dt, 0.1)
      
      res <- residuals(fit)
      mean_x <- mean(xlag_s, na.rm = TRUE)
      sigma_hat <- sqrt(var(res) / (max(mean_x, 1e-6) * dt))
      sigma_j[s] <- max(sigma_hat, 0.01)
    }
    
    # --------- build Gamma based on ORIGINAL labels ----------
    Gamma <- matrix(NA_real_, states, states)
    for (i in 1:states) {
      idx_i <- which(cluster[1:(TT - 1)] == i)
      if (length(idx_i) == 0) {
        Gamma[i, ] <- 1 / states
        next
      }
      num <- numeric(states)
      for (k in 1:states) {
        N_ik <- sum(cluster[idx_i + 1] == k)
        num[k] <- max(N_ik, 1e-2)
      }
      Gamma[i, ] <- num / sum(num)
    }
    
    # ---------- ordering rule (states = 2 only; calm=lower theta by default) ----------
    if (states == 2) {
      # 你现在逻辑：theta_1 > theta_2 -> swap
      # 这等价于让 theta 小的当 state 1（更 calm）
      if (theta_j[1] > theta_j[2]) {
        cluster <- 3 - cluster
        kappa_j <- rev(kappa_j)
        theta_j <- rev(theta_j)
        sigma_j <- rev(sigma_j)
        Gamma <- Gamma[2:1, 2:1]
      }
    }
    
    # ---- store ----
    cluster_mat[j, ] <- cluster
    kappa_mat[, j] <- kappa_j
    theta_mat[, j] <- theta_j
    sigma_mat[, j] <- sigma_j
    Gamma_list[[j]] <- Gamma
  }
  
  # =========================
  # Build param_mat: (states*d) x 3, regime1 x1..xd, regime2 x1..xd, ...
  # =========================
  param_mat <- cbind(
    kappa = as.vector(t(kappa_mat)),
    theta = as.vector(t(theta_mat)),
    sigma = as.vector(t(sigma_mat))
  )
  colnames(param_mat) <- c("kappa", "theta", "sigma")
  
  # ---- estimate R using standardized innovations ----
  R_hat <- estimate_R_sample(
    X = X,
    Reg = cluster_mat,     # use inferred regimes
    param_mat = param_mat,
    dt = dt
  )
  
  list(
    cluster    = cluster_mat,
    reg_param  = param_mat,
    Gamma_list = Gamma_list,
    R_hat      = R_hat
  )
}




init <- initial_heuristic_multi(X, states = 2, dt = 1/252)
plot_multi_series(init$cluster)
plot_multi_series(Reg)
init$reg_param
init$Gamma_list
init$R_hat
R
# initial_heuristic(X[1,],2, T)

# init$cluster













estimate_R_sample <- function(X, Reg = NULL, param_mat, dt = 1/252, eps = 1e-8) {

  
  d <- nrow(X)
  TT <- ncol(X) - 1  # number of increments
  
  # increments
  dX <- X[, 2:(TT+1), drop=FALSE] - X[, 1:TT, drop=FALSE]
  Xlag <- X[, 1:TT, drop=FALSE]
  
  # helper to get (kappa, theta, sigma) for each coordinate/time
  get_pars_at_k <- function(j, k) {
    if (is.null(Reg)) {
      # no regime provided: use regime 1 as default
      row <- j
    } else {
      # Reg could be d x TT or d x (TT+1); use k
      s_jk <- Reg[j, min(k, ncol(Reg)), drop=TRUE]  # 0/1 or 1/2 depends on your sim
      # if your sim uses 0/1, map to 1/2
      if (s_jk %in% c(0, 1)) s_jk <- s_jk + 1
      row <- (s_jk - 1) * d + j
    }
    list(
      kappa = param_mat[row, "kappa"],
      theta = param_mat[row, "theta"],
      sigma = param_mat[row, "sigma"]
    )
  }
  
  # build standardized shocks eps_{j,k} ≈ (ΔX - kappa(theta-X)dt) / (sigma*sqrt(X)*sqrt(dt))
  E <- matrix(NA_real_, nrow = TT, ncol = d)  # TT x d
  for (k in 1:TT) {
    for (j in 1:d) {
      pars <- get_pars_at_k(j, k)
      mu_jk <- pars$kappa * (pars$theta - Xlag[j, k]) * dt
      denom <- pars$sigma * sqrt(pmax(Xlag[j, k], eps)) * sqrt(dt)
      E[k, j] <- (dX[j, k] - mu_jk) / denom
    }
  }
  
  # sample correlation of standardized shocks
  R_hat <- cor(E, use = "pairwise.complete.obs")
  
  # (optional) small numerical fix to ensure symmetry and ones on diagonal
  R_hat <- (R_hat + t(R_hat)) / 2
  diag(R_hat) <- 1
  
  R_hat
}


dt <- 1/252
R_hat <- estimate_R_sample(X = X, Reg = Reg, param_mat = param_mat, dt = dt)
print(R_hat)



initial_heuristic(X[1,], 2)
