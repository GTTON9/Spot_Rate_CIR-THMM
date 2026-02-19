calculate_OM_for_likelihood <- function(result_with_paths,
                                        kappa, theta, sigma,
                                        H = 5/252, N = 100,
                                        lambda_om = 1.0,
                                        use_upper = TRUE,
                                        use_lower = TRUE) {
  T_periods <- nrow(result_with_paths)
  nstates <- length(kappa)
  
  log_b <- matrix(0, nrow = nstates, ncol = T_periods)
  
  for (t in 1:T_periods) {
    a_upper <- result_with_paths$path_upper_a[t]
    b_upper <- result_with_paths$path_upper_b[t]
    a_lower <- result_with_paths$path_lower_a[t]
    b_lower <- result_with_paths$path_lower_b[t]
    
    for (i in 1:nstates) {
      Li_t <- 0
      
      if (use_upper && !is.na(a_upper) && !is.na(b_upper)) {
        Li_t <- Li_t + calculate_OM_single(
          a = a_upper, b = b_upper,
          kappa = kappa[i], theta = theta[i], sigma = sigma[i],
          H = H, N = N
        )
      }
      
      if (use_lower && !is.na(a_lower) && !is.na(b_lower)) {
        Li_t <- Li_t + calculate_OM_single(
          a = a_lower, b = b_lower,
          kappa = kappa[i], theta = theta[i], sigma = sigma[i],
          H = H, N = N
        )
      }
      
      # paper: log b_i(Ot) = - L_i[phi] (up to additive const independent of state)
      log_b[i, t] <- -lambda_om * Li_t
    }
  }
  
  log_b
}