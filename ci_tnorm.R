
g_log <- function() {
  list(
    g     = function(x) log(x),
    g_inv = function(x) exp(x),
    g_der = function(x) 1 / x
  )
}

g_exp <- function() {
  list(
    g     = function(x) exp(x),
    g_inv = function(x) log(x),
    g_der = function(x) exp(x)
  )
}

g_identity <- function() {
  list(
    g     = function(x) x,
    g_inv = function(x) x,
    g_der = function(x) 1
  )
}

# logit transform: g(x) = log(x / (1 - x)), for x in (0, 1)
# use for probability parameter
g_logit <- function() {
  list(
    g     = function(x) log(x / (1 - x)),
    g_inv = function(x) 1 / (1 + exp(-x)),          # inverse logit / expit
    g_der = function(x) 1 / (x * (1 - x))
  )
}

# arctangent transform: g(x) = atan(x)
g_arctan <- function() {
  list(
    g     = function(x) atan(x),
    g_inv = function(x) tan(x),
    g_der = function(x) 1 / (1 + x^2)
  )
}

ci_TNORM_Jeng <- function(t, t_c, g, alpha = 0.1) {
  # g: a list with elements g, g_inv, g_der (e.g. g_log(), g_exp(), g_identity())
  
  # maximum likelihood estimations
  fit <- mle_sev(t, t_c)
  
  mle <- fit$par
  mu_hat <- mle[1]
  sigma_hat <- mle[2]
  
  cov_matrix <- fit$fisher
  se <- sqrt(diag(cov_matrix))
  se_mu <- se[1]
  se_sigma <- se[2]
  
  z <- qnorm(p = 1 - alpha / 2)
  
  # build a TNORM CI for a scalar estimate theta_hat with se se_theta,
  # using transformation g: CI on g-scale is g(theta_hat) +- z*se(g(theta_hat)),
  # then back-transformed via g_inv
  tnorm_ci <- function(theta_hat, se_theta) {
    g_hat <- g$g(theta_hat)
    se_g  <- g$g_der(theta_hat) * se_theta
    ci_g  <- c(g_hat - z * se_g, g_hat + z * se_g)
    g$g_inv(ci_g)
  }
  
  ci_mu    <- tnorm_ci(mu_hat, se_mu)
  ci_sigma <- tnorm_ci(sigma_hat, se_sigma)
  
  # Delta method standard error (on the original t_p scale)
  delta_se <- function(p) {
    grad <- c(1, log(-log(1 - p)))
    sqrt(t(grad) %*% cov_matrix %*% grad)
  }
  
  # t_0.1
  t01_hat <- mu_hat + sigma_hat * log(-log(1 - 0.1))
  se_t01  <- delta_se(0.1)
  ci_t01  <- tnorm_ci(t01_hat, se_t01)
  
  # t_0.5
  t05_hat <- mu_hat + sigma_hat * log(-log(1 - 0.5))
  se_t05  <- delta_se(0.5)
  ci_t05  <- tnorm_ci(t05_hat, se_t05)
  
  return(list(
    mu    = unname(ci_mu),
    sigma = unname(ci_sigma),
    t_0.1 = unname(ci_t01),
    t_0.5 = unname(ci_t05)
  ))
}

ci_TNORM_Jeng(t = ball_bearing_data, t_c = 60, g = g_arctan(), alpha = 0.1)[1]
ci_NORM_Jeng(t = ball_bearing_data, t_c = 60, alpha = 0.1)[1]
