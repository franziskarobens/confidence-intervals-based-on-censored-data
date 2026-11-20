library(plotly)


mu_vals <- seq(-10, 10, length.out = 50)
sigma_vals <- seq(0.1, 2, length.out = 50)

loglik_matrix <- matrix(NA, 
                        nrow = length(mu_vals), 
                        ncol = length(sigma_vals))


for(i in 1:length(mu_vals)) {
  for(j in 1:length(sigma_vals)) {
    loglik_matrix[i, j] <- loglik_sev(
      par = c(mu_vals[i], sigma_vals[j]),
      t = t,
      delta = delta,
      t_c = t_c
    )
  }
}


# heatmap
image(mu_vals,
      sigma_vals,
      loglik_matrix,
      xlab = "mu",
      ylab = "sigma")


