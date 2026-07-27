
ci_NORM_Jeng <- function(t, t_c, alpha = 0.1) {
  # Jeng NORM, 
  
  # maximum likelihood estimations
  fit <- mle_sev(t, t_c)
  
  mle <- fit$par
  mu_hat <- mle[1]
  sigma_hat <- mle[2]
  
  cov_matrix <- fit$fisher
  se <- sqrt(diag(cov_matrix))
  se_mu <- se[1]
  se_sigma <- se[2]
  # print(cov_matrix)
  # cat("se mu = ", se_mu, " se_sigma = ", se_sigma, "\n")
  
  z <- qnorm(p = 1 - alpha/2)
  
  ci_mu <- c(mu_hat - z * se_mu, mu_hat + z * se_mu)
  ci_sigma <- c(sigma_hat - z * se_sigma, sigma_hat + z * se_sigma)
  
  # Delta method standard error
  delta_se <- function(p) {
    g <- c(1, log(-log(1 - p)))
    sqrt(t(g) %*% cov_matrix %*% g)
  }
  
  # TODO mit den t_p CI stimmt noch etwas nicht
  # t_0.1
  t01_hat <- mu_hat + sigma_hat * log(- log(1 - 0.1))
  se_t01 <- delta_se(0.1)
  ci_t01 <- c(t01_hat - z * se_t01, t01_hat + z * se_t01)
  
  # t_0.5
  t05_hat <- mu_hat + sigma_hat * log(- log(1 - 0.5))
  se_t05 <- delta_se(0.5)
  ci_t05 <- c(t05_hat - z * se_t05, t05_hat + z * se_t05)
  
  return(list(
    mu = unname(ci_mu),
    sigma = unname(ci_sigma),
    t_0.1 = ci_t01,
    t_0.5 = ci_t05
  ))
}






mle_Sundberg <- function(t, t_c) {
  N <- sum(1 - get_delta(data = t, t_c = t_c))
  
  return(sum(t) / N)
}


# TODO warum mit dem mle von sundberg komplett anderes ergebnis als mit dem 
# von sev_mle ???
ci_NORM_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  # maximum likelihood estimations
  mu_hat <- mle_Sundberg(t, t_c)
  
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  z <- qnorm(p = 1 - alpha / 2)
  
  ci <- c(mu_hat * (1 - z / sqrt(N)), mu_hat * (1 + z / sqrt(N)))
  
  return(unname(ci))
}


# TODO warum sind die CI so unterschiedlich? Liegt es an der Unterschiedlichen 
# zugrundeliegenden verteilung (SEV)?
ci_NORM_Jeng(t = ball_bearing_data, t_c = 40, alpha = 0.1)[1]
ci_NORM_Jeng(t = ball_bearing_data, t_c = 60, alpha = 0.1)[1]

ci_NORM_Sundberg(t = ball_bearing_data, t_c = 40, alpha = 0.1)  
ci_NORM_Sundberg(t = ball_bearing_data, t_c = 60, alpha = 0.1)  



