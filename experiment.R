


N <- 252
T <- 1

S0 <- 100
v0 <- 0.1

n_days <- 252
n_intraday <- 2400
nstates <- 2


n_dim <- 4
Gamma <- matrix(
  c(0.99, 0.01,
    0.01, 0.99),
  nrow = 2
)

mu <- c(0.3, 0.3)
kappa <- c(10, 5)
theta <- c(0.1, 0.5)
sigma <- c(0.05, 0.05)
rho <- c(-0.9, -0.9)

regime_params <- list(
  list(mu = mu[1], kappa = kappa[1], theta = theta[1], sigma = sigma[1], rho = rho[1]),
  list(mu = mu[2], kappa = kappa[2], theta = theta[2], sigma = sigma[2], rho = rho[2])
)


Reg_chain <-simulate_Regime(Gamma, N, n_dim)


simulted_Heston_Series <- simulate_heston(
  S0, v0,
  Reg_chain,
  regime_params,
  T = 1,
  N = N,
  M = 1,
  method = "E",
  interp = TRUE,
  substeps = n_intraday
)
n_intraday <- 2400

V_daily <- simulted_Heston_Series$V_paths[, seq(1, ncol(simulted_Heston_Series$V_paths), by = n_intraday)]
S_daily <- simulted_Heston_Series$S_paths[, seq(1, ncol(simulted_Heston_Series$S_paths), by = n_intraday)]
plot_multi_series(Reg_chain)
plot_multi_series(V_daily)
plot_multi_series(S_daily)





