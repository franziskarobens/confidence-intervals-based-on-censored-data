# Monte Carlo coverage probability study

number_simulations <- 2000


Er_vals <- c(3, 5, 7, 10, 15, 20, 30, 50, 100)
pf_vals <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

alpha <- c(0.025, 0.05)
alpha <- 0.05



# ? se <- sqrt(alpha * (1 - alpha) / monte_carlo_samples)

# initialize matrix with coverage probabilities
coverage_probabilities <- matrix(0, 
                                 nrow = length(Er_vals), ncol = length(pf_vals))
rownames(coverage_probabilities) <- Er_vals
colnames(coverage_probabilities) <- pf_vals

r_count <- matrix(0, nrow = length(Er_vals), ncol = length(pf_vals))
rownames(r_count) <- Er_vals
colnames(r_count) <- pf_vals


for (i in seq_along(Er_vals)) {
  for (j in seq_along(pf_vals)) {
    
    Er <- Er_vals[i]
    pf <- pf_vals[j]
    
    n <- ceiling(Er / pf)
    
    t_c <- qsev(q = pf)
    
    is_in_count <- 0
    
    for (sim in 1:number_simulations) {
      
      # generate dataset
      uncensored_sample <- rsev(n = n)
      
      # censoring
      censored_sample <- pmin(uncensored_sample, t_c)
      delta <- get_delta(uncensored_sample, t_c = t_c)
      
      # calculate number of failures before t_c
      r <- sum(delta)
      
      # get confidence interval 
      ci <- ci_NORM(t = uncensored_sample, t_c = t_c, alpha = alpha)
      
      # check if parameter is in CI
      if (ci[1] <= real_parameter && real_parameter <= ci[2]) {
        is_in_count <- is_in + 1
      }
      
      cat("sim=", sim, " r=", r, " ci=", ci, " is_in=", is_in, "\n")
    }
    
    coverage_probabilities[i, j] <- is_in / number_simulations
  }
}



