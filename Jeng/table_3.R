p_f_vals <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9)
E_r_vals <- c(3, 5, 7, 10, 15, 20, 30, 50, 100)


##### Berechne die theoretischen Werte in den Klammern von Table 3

p_f <- 0.1
E_r <- 7

n <- E_r / p_f

prob_smaller_1 <- (1 - p_f)^(n-1) * (1 - p_f + n * p_f)
2000 * prob_smaller_1 # value in braces


#### Simulate with Monte Carlo

M <- 2000

tab_sim <- sapply(p_f_vals, function(p) {
  sapply(E_r_vals, function(E) {
    n <- round(E / p)
    r <- rbinom(M, size = n, prob = p)     # Anzahl Ausfälle pro Simulation
    sum(r <= 1)
  })
})
