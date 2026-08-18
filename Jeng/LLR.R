
ci_LLR <- function(t, t_c, alpha) {
  
  # maximum likelihood estimations
  fit <- mle_sev(t, t_c)
  
  mu_hat <- fit$par[1]
  sigma_hat <- fit$par[2]
  
  cov_matrix <- fit$fisher
  se <- sqrt(diag(cov_matrix))
  se_mu <- se[1]
  se_sigma <- se[2]
  
  
  loglik_max <- loglik_sev(fit$par, t, get_delta(t, t_c), t_c)
  
  cutoff <- loglik_max - 0.5 * qchisq(1 - alpha, df = 1)
  
  # grid for sigma
  sigma_grid <- seq(0.1 * sigma_hat, 3 * sigma_hat, length = 500)
  
  loglik_sigma <- sapply(sigma_grid, function(s) {
    
    opt <- optim(
      par = mu_hat,
      fn = function(mu) -loglik_sev(c(mu, s), t, t_c, delta = get_delta(data = t, t_c = t_c))
    )
    
    loglik_sev(c(opt$par, s), t, t_c, delta = get_delta(data=t, t_c = t_c))
  })
  
  sigma_ci <- range(sigma_grid[loglik_sigma >= cutoff])
  
  # quantile function
  tp <- function(p, mu, sigma) {
    mu + sigma * log(-log(1 - p))
  }
  
  # t_0.1 CI
  t01_vals <- sapply(1:length(sigma_grid), function(i) {
    
    s <- sigma_grid[i]
    
    opt <- optim(
      par = mu_hat,
      fn = function(mu) - loglik_sev(c(mu, s), t, t_c, delta = get_delta(data = t, t_c = t_c))
    )
    
    tp(0.1, opt$par, s)
  })
  
  t01_ci <- range(t01_vals[loglik_sigma >= cutoff])
  
  # t_0.5 CI
  t05_vals <- sapply(1:length(sigma_grid), function(i) {
    
    s <- sigma_grid[i]
    
    opt <- optim(
      par = mu_hat,
      fn = function(mu) - loglik_sev(c(mu, s), t, t_c, delta = get_delta(data=t, t_c = t_c))
    )
    
    tp(0.5, opt$par, s)
  })
  
  t05_ci <- range(t05_vals[loglik_sigma >= cutoff])
  
  return(list(
    sigma = sigma_ci,
    t_0.1 = t01_ci,
    t_0.5 = t05_ci
  ))
}

ci_LLR(t = ball_bearing_data, t_c = 60, alpha = 0.9)

