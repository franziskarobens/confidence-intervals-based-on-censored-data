test_get_delta <- function() {
  
  data <- c(1, 2, 3, 4, 5, 6, 2, 1, 4, 7, 9)
  t_c <- 5
  
  delta_expected <- c(1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0)
  delta_actual <- get_delta(data, t_c)
  
  
  if (!isTRUE(all.equal(delta_expected, delta_actual, tolerance = 1e-8))) {
    stop("\nExpected delta = ", 
         paste(delta_expected, collapse = ", "), 
         "\nGot delta = ", 
         paste(delta_actual, collapse = ", "))
  } else {
    cat("\033[32mTest passed for get_delta()\033[39m\n")
  }
  
}


test_lik_sev <- function(){
  
  par <- c(1, 2)
  mu <- par[1]
  sigma <- par[2]
  
  t <- c(1, 4, 6)
  t_c <- 5
  delta <- c(1, 1, 0)
  
  lik_expected <- dsev(log(1), mu = mu, sigma = sigma) * dsev(log(4), mu = mu, sigma = sigma) * (1 - psev(log(5), mu = mu, sigma = sigma))
  lik_actual <- lik_sev(par = par, t = t, delta = delta, t_c = t_c)
  
  if (!isTRUE(all.equal(lik_expected, lik_actual, tolerance = 1e-8))) {
    stop("\nExpected likelihood = ", 
         paste(lik_expected, collapse = ", "), 
         "\nGot likelihood = ", 
         paste(lik_actual, collapse = ", "))
  } else {
    cat("\033[32mTest passed for lik_sev()\033[39m\n")
  }
  
}

test_loglik_sev <- function() {
  t <- c(1, 4, 6)
  delta <- c(1, 1, 0)
  t_c <- 5
  
  par <- c(1, 2)
  mu <- par[1]
  sigma <- par[2]
  
  loglik_expected <- log(dsev(log(1), mu = mu, sigma = sigma)) + log(dsev(log(4), mu = mu, sigma = sigma)) + log(1 - psev(log(5), mu = mu, sigma = sigma))
  loglik_actual <- loglik_sev(par = par, t = t, delta = delta, t_c = t_c)
  
  if (!isTRUE(all.equal(loglik_expected, loglik_actual, tolerance = 1e-8))) {
    stop("\nExpected log likelihood = ", 
         paste(loglik_expected, collapse = ", "), 
         "\nGot log likelihood = ", 
         paste(loglik_actual, collapse = ", "))
  } else {
    cat("\033[32mTest passed for loglik_sev()\033[39m\n")
  }
}


test_sev_functions <- function() {
  y <- 3.5
  mu = 4
  sigma = 9
  
  # dsev
  temp <- (y - mu) / sigma
  expected_dsev <- (1/sigma) * exp(temp - exp(temp))
  actual_dsev <- dsev(y, mu, sigma)
  
  # psev
  expected_psev <- 1 - exp(- exp(temp))
  actual_psev <- psev(y, mu, sigma)
  
  # check dsev
  if (!isTRUE(all.equal(expected_dsev, actual_dsev, tolerance = 1e-8))) {
    stop("\nExpected dsev = ", 
         paste(expected_dsev, collapse = ", "), 
         "\nGot dsev = ", 
         paste(actual_dsev, collapse = ", "))
  } else {
    cat("\033[32mTest passed for dsev()\033[39m\n")
  }
  
  # check psev
  if (!isTRUE(all.equal(expected_psev, actual_psev, tolerance = 1e-8))) {
    stop("\nExpected psev = ", 
         paste(expected_psev, collapse = ", "), 
         "\nGot psev = ", 
         paste(actual_psev, collapse = ", "))
  } else {
    cat("\033[32mTest passed for psev()\033[39m\n")
  }
  
}


test_mle_sev <- function() {
  ## TODO 
}

test_sev_functions()
test_get_delta()
test_lik_sev()
test_loglik_sev()
