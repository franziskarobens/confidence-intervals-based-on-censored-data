# SEV density
dsev <- function(y, mu = 0, sigma = 1){
  z <- (y - mu)/sigma
  1/sigma * exp(z - exp(z))
}

# SEV CDF
psev <- function(y, mu = 0, sigma = 1){
  z <- (y - mu)/sigma
  1 - exp(-exp(z))
}

# random samples
rsev <- function(n, mu = 0, sigma = 1) {
  # Inverse Transform Sampling
  u <- runif(n)
  mu + sigma * log(-log(1 - u))
}

# quantiles from sev dist'n
qsev <- function(q, mu = 0, sigma = 1) {
  c_q <- log(- log(1 - q))
  mu + sigma * c_q 
}



plot_sev_density <- function(mu, sigma) {
  # Create range of y values
  y <- seq(..., ..., length.out = 500)
  
  # Compute density
  density <- dsev(y, mu = mu, sigma = sigma)
  
  # Plot
  plot(y, density, 
       type = "l",
       lwd = 2,
       xlab = "y",
       ylab = "Density",
       main = "SEV Density (mu = ..., sigma = ...)")
}




lik_sev <- function(par, t, delta, t_c) {
  
  mu <- par[1]
  sigma <- par[2]
  
  t_obs <- pmin(t, t_c)
  log_t <- log(t_obs)
  
  uncensored <- dsev(y = log_t, mu = mu, sigma = sigma)
  censored <- 1 - psev(y = log_t, mu = mu, sigma = sigma)
  all_factors <- uncensored^delta * censored^(1 - delta)
  return(prod(all_factors))
}



loglik_sev <- function(par, t, delta, t_c) {
  
  mu <- par[1]
  sigma <- par[2]
  
  r <- sum(1 - delta) # number of censored data points
  
  uncensored <- delta * log(dsev(y = log(t), mu = mu, sigma = sigma))
  censored <- r * log(1 - psev(y = log(t_c), mu = mu, sigma = sigma))
  
  loglik <- sum(uncensored) + censored
  
  # Handle NaN / Inf
  if (!is.finite(loglik)) {
    loglik <- -1e+10   # very bad likelihood
  }
  
  # cat("loglik =", loglik, "mu =", mu, "sigma =", sigma, "\n")
  return(loglik)
}



mle_sev <- function(t, t_c, draw = FALSE) {
  
  delta <- get_delta(data = t, t_c = t_c)
  
  if(draw) {
    
    # ---- Fix sigma = 1, vary mu ----
    mu_vals <- seq(mean(y) - 2*sd(t), mean(t) + 2*sd(t), length.out = 100)
    
    ll_mu <- sapply(mu_vals, function(mu){
      loglik_sev(par = c(mu, 1), t = t, delta = delta, t_c = t_c)
    })
    
    plot(mu_vals, ll_mu,
         type = "l",
         lwd = 2,
         xlab = expression(mu),
         ylab = "Log-Likelihood",
         main = expression("Log-Likelihood vs " * mu * " (sigma = 1)"))
    
    
    # ---- Fix mu = 0, vary sigma ----
    sigma_vals <- seq(0.1, 3*sd(t), length.out = 100)
    
    ll_sigma <- sapply(sigma_vals, function(sigma){
      loglik_sev(par = c(0, sigma), t = t, delta = delta, t_c = t_c)
    })
    
    plot(sigma_vals, ll_sigma,
         type = "l",
         lwd = 2,
         xlab = expression(sigma),
         ylab = "Log-Likelihood",
         main = expression("Log-Likelihood vs " * sigma * " (mu = 0)"))
  }
  
  
  ## TODO maximizing doenst fully work, gets stuck in infinity
  fit <- optim(
    par = c(mu = mean(t), sigma = sd(t)), # starting values
    fn = function(...) -loglik_sev(...), # optim() minimizes
    t = t,
    t_c = t_c,
    delta = delta,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-6), 
    hessian = TRUE
  )
  
  if(fit$convergence != 0) return(NA)
  
  fit$fisher <- solve(fit$hessian)
  
  return(fit)
}

