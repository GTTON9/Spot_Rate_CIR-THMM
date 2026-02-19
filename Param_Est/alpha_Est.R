as_design <- function(X, r) { # create the design matrix
  X <- as.matrix(X)
  r <- as.numeric(r)
  
  if (nrow(X) == 3 && length(r) == ncol(X)) {
    Xn <- t(X)
  } else if (ncol(X) == 3 && length(r) == nrow(X)) {
    Xn <- X
  } else {
    stop("Dimension mismatch: expect X to be 3xN (or Nx3) and r length N.")
  }
  colnames(Xn) <- paste0("X", 1:ncol(Xn))
  list(X = Xn, r = r)
}

estimate_alpha_ols <- function(X, r) {
  dat <- as_design(X, r)
  fit <- lm(dat$r ~ dat$X)
  coef(fit)  
}



alpha_est <- estimate_alpha_ols(X, r)

