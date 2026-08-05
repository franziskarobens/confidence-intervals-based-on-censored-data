
# construct NORM confidence interval regarding to Jeng (3.1) for SEV distribution
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








ci_NORM_Jeng(t = ball_bearing_data, t_c = 40, alpha = 0.1)[1]
ci_NORM_Jeng(t = ball_bearing_data, t_c = 60, alpha = 0.1)[1]

C1_Sundberg(t = ball_bearing_data, t_c = 40, alpha = 0.1)  
C1_Sundberg(t = ball_bearing_data, t_c = 60, alpha = 0.1) 


# the CI differ a lot because the assumption of the underlying dist'n is different 
# (SEV for Jeng vs. Exp for Sundberg)