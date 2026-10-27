
# construct NORM confidence interval regarding to Jeng (3.1) for SEV distribution
NORM_CI <- function(t, t_c, alpha = 0.1) {
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
  
  
  # t_0.1
  c_01 <- log(- log(1 - 0.1))
  hat_y_01 <- mu_hat + c_01 * sigma_hat
  hat_t_01 <- exp(hat_y_01)
  se_t01 <- delta_se(0.1)
  ci_t01 <- c(hat_t_01 - z * se_t01, hat_t_01 + z * se_t01)

  # t_0.5
  c_05 <- log(- log(1 - 0.5))
  hat_y_05 <- mu_hat + c_05 * sigma_hat
  hat_t_05 <- exp(hat_y_05)
  se_t05 <- delta_se(0.5)
  ci_t05 <- c(hat_t_05 - z * se_t05, hat_t_05 + z * se_t05)
  
  return(list(
    mu = unname(ci_mu),
    sigma = unname(ci_sigma),
    t_0.1 = ci_t01,
    t_0.5 = ci_t05
  ))
}






NORM_CI(t = ball_bearing_data, t_c = 40, alpha = 0.1)[1]
NORM_CI(t = ball_bearing_data, t_c = 60, alpha = 0.1)[1]

r_40 <- sum(ball_bearing_data <= 40)
MLE_40 <- mle_sev(t = ball_bearing_data, t_c = 40)$par["mu"]
C1 <- C1_Sundberg(r = r_40, alpha = 0.1) * MLE_40
C1

r_60 <- sum(ball_bearing_data <= 60)
MLE_60 <- mle_sev(t = ball_bearing_data, t_c = 60)$par["mu"]
C2 <- C1_Sundberg(r = r_60, alpha = 0.1) * MLE_60
C2

# the CI differ a lot because the assumption of the underlying dist'n is different 
# (SEV for Jeng vs. Exp for Sundberg)


# Table 2:

NORM_CI(t = ball_bearing_data, t_c = 40, alpha = 0.1)
# Paper: [0.00, 0.89], [17.97, 52.73], [7.69, 148.92]

NORM_CI(t = ball_bearing_data, t_c = 60, alpha = 0.1)
# Paper: [0.15, 0.50], [2.146, 44.73], [49.02, 72.85]
