plot_multi_series <- function(Reg_mat,
                               xlab = "Time",
                               ylab = "Value") {
  
  stopifnot(is.matrix(Reg_mat))
  
  n_dim <- nrow(Reg_mat)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(n_dim, 1), mar = c(3,4,2,1))
  
  for (j in 1:n_dim) {
    plot(
      Reg_mat[j, ],
      type = "l",               
      xlab = xlab,
      ylab = ylab
    )
    grid()
  }
}
