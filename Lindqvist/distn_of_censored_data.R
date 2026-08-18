# Description ----------------------------------------------

# Simulates with Monte Carlo the asymptotic distribution of 
# standardized MLEs for the exponential rate parameter λ under different 
# censoring schemes (uncensored, Type I, Type II, Hybrid Type I, Hybrid Type II)

# Helper Function ----------------------------------------------

plot_standardized_MLE <- function(MLEs, lambda, n_s, C = NULL, is_censored) {
  
  if (is_censored) {
    title <- paste0("Standardized censored MLE vs. N(0,1) using t ~ Exp(", lambda, ") and C = ", C)
    expression_xlab <- expression(sqrt(N) * (hat(lambda) - lambda) / hat(lambda))
    if (is.null(C)) {
      stop("If data is censored, give C as an argument to plot_standardized_MLE()")
    }
  } else {
    title <- paste0("Standardized uncensored MLE vs. N(0,1) using t ~ Exp(", lambda, ")")
    expression_xlab <- expression(sqrt(n) * (hat(lambda) - lambda) / hat(lambda))
  }
  
  standardized <- sqrt(n_s) * (MLEs - lambda) / MLEs

  # Histogram (auf Dichteskala, damit die Kurven vergleichbar sind)
  hist(standardized, breaks = 50, freq = FALSE,
      main = title,
      xlab = expression_xlab
    )

  # Normalverteilung (rot) drüberlegen
  curve(dnorm(x, mean = 0, sd = 1), col = "red", lwd = 2, add = TRUE)

  # Kerndichteschätzung (blau) drüberlegen
  lines(density(standardized), col = "blue", lwd = 2)

  # small p value -> hint for non N(0, 1) distribution
  p_value <- ks.test(standardized, "pnorm", mean = 0, sd = 1)$p.value

  legend("topright",
         legend = c("N(0,1)", "Kernel density",
                    paste0("KS p-value = ", signif(p_value, 3))),
         col = c("red", "blue", NA),
         lwd = c(2, 2, NA),
         bty = "n")
}


# Type I Censoring ----------------------------------------------

lambda <- 10
C <- 0.1
n <- 1000
n_simulations <- 100000

run_censoring_simulation_type_I <- function(n_simulations, n, C) {
  MLEs_u <- numeric(n_simulations)
  MLEs_c <- numeric(n_simulations)
  N_u <- rep(n, n_simulations)
  N_c <- numeric(n_simulations)

  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  on.exit(close(pb))  # stellt sicher, dass die Progress Bar auch bei einem Fehler geschlossen wird

  for (i in 1:n_simulations) {
    t <- rexp(n = n, rate = lambda)
    t_cens <- pmin(t, C)
    N <- sum(t <= C) # number of uncensored observations

    MLEs_u[i] <- n / sum(t)
    MLEs_c[i] <- N / sum(t_cens)
    N_c[i] <- N

    setTxtProgressBar(pb, i)
  }

  list(
    uncensored_MLEs = MLEs_u, censored_MLEs = MLEs_c, 
    uncensored_N = N_u, censored_N = N_c
  )
}

res_exp <- run_censoring_simulation_type_I(
  n_simulations = n_simulations,
  n = n,
  C = C
)

plot_standardized_MLE(MLEs = res_exp$uncensored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$uncensored_N,
  is_censored = FALSE
)

plot_standardized_MLE(MLEs = res_exp$censored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$censored_N,
  C = C,
  is_censored = TRUE
)


# Type II Censoring ----------------------------------------------

lambda <- 10

m <- 650

n <- 1000
n_simulations <- 10000

run_censoring_simulation_type_II <- function(n_simulations, n, C, m) {
  if (m >= n) {
    stop("m must be smaller than n. m = ", m, " and n = ", n)
  }


  MLEs_u <- numeric(n_simulations)
  MLEs_c <- numeric(n_simulations)
  N_u <- rep(n, n_simulations)
  N_c <- rep(m, n_simulations)

  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  on.exit(close(pb))

  for (i in 1:n_simulations) {
    t <- rexp(n, rate = lambda)

    t_sorted <- sort(t)
    t_m <- t_sorted[m]

    t_cens <- c(t_sorted[1:m], rep(t_m, n - m))

    MLEs_u[i] <- n / sum(t)
    MLEs_c[i] <- m / sum(t_cens)

    setTxtProgressBar(pb, i)
  }

  list(
    uncensored_MLEs = MLEs_u, censored_MLEs = MLEs_c, 
    uncensored_N = N_u, censored_N = N_c
  )
}

res_exp <- run_censoring_simulation_type_II(
  n_simulations = n_simulations,
  n = n,
  C = C, 
  m = m
)

plot_standardized_MLE(MLEs = res_exp$uncensored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$uncensored_N,
  is_censored = FALSE
)

plot_standardized_MLE(MLEs = res_exp$censored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$censored_N,
  C = C,
  is_censored = TRUE
)


# Hybrid Type I Censoring ----------------------------------------------

lambda <- 10

C <- 0.1
m <- 650

n <- 1000
n_simulations <- 10000

run_censoring_simulation_hybrid_type_I <- function(n_simulations, n, C, m) {
  if (m >= n) {
    stop("m must be smaller than n. m = ", m, " and n = ", n)
  }


  MLEs_u <- numeric(n_simulations)
  MLEs_c <- numeric(n_simulations)
  N_u <- rep(n, n_simulations)
  N_c <- numeric(n_simulations)
  censoring_case <- numeric(n_simulations)

  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  on.exit(close(pb))

  for (i in 1:n_simulations) {
    t <- rexp(n, rate = lambda)

    if (all(t < C)) {
      warning("No time censoring occured. Parameter C might be too high. C = ", C)
    }

    t_sorted <- sort(t)
    t_m <- t_sorted[m]

    if (t_m < C) { # exactly m observations
      t_cens <- c(t_sorted[1:m], rep(t_m, n - m))
      N <- m
      censoring_case[i] <- 0
    } else {
      t_cens <- pmin(t, C)
      N <- sum(t <= C)
      censoring_case[i] <- 1
    }

    MLEs_u[i] <- n / sum(t)
    MLEs_c[i] <- N / sum(t_cens)
    N_c[i] <- N

    setTxtProgressBar(pb, i)
  }

  list(
    uncensored_MLEs = MLEs_u, censored_MLEs = MLEs_c, 
    uncensored_N = N_u, censored_N = N_c, 
    censoring_case = censoring_case
  )
}

res_exp <- run_censoring_simulation_hybrid_type_I(
  n_simulations = n_simulations,
  n = n,
  C = C, 
  m = m
)

plot_standardized_MLE(MLEs = res_exp$uncensored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$uncensored_N,
  is_censored = FALSE
)

plot_standardized_MLE(MLEs = res_exp$censored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$censored_N,
  C = C,
  is_censored = TRUE
)

hist(res_exp$censored_N, breaks = 30)
# gets cut off at m, otherwise normal distributed around WHAT VALUE?



# Hybrid Type II Censoring ----------------------------------------------

lambda <- 10

C <- 0.1
m <- 650

n <- 1000
n_simulations <- 10000

run_censoring_simulation_hybrid_type_II <- function(n_simulations, n, C, m) {
  if (m >= n) {
    stop("m must be smaller than n. m = ", m, " and n = ", n)
  }


  MLEs_u <- numeric(n_simulations)
  MLEs_c <- numeric(n_simulations)
  N_u <- rep(n, n_simulations)
  N_c <- numeric(n_simulations)

  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  on.exit(close(pb))

  for (i in 1:n_simulations) {
    t <- rexp(n, rate = lambda)

    if (all(t < C)) {
      warning("No time censoring occured. Parameter C might be too high. C = ", C)
    }

    t_sorted <- sort(t)
    t_m <- t_sorted[m]

    if (t_m > C) { # exactly m observations
      t_cens <- c(t_sorted[1:m], rep(t_m, n - m))
      N <- m
    } else { # first m observations and then time stop
      t_cens <- pmin(t, C)
      N <- sum(t <= C)
    }

    MLEs_u[i] <- n / sum(t)
    MLEs_c[i] <- N / sum(t_cens)
    N_c[i] <- N

    setTxtProgressBar(pb, i)
  }

  list(
    uncensored_MLEs = MLEs_u, censored_MLEs = MLEs_c, 
    uncensored_N = N_u, censored_N = N_c
  )
}

res_exp <- run_censoring_simulation_hybrid_type_II(
  n_simulations = n_simulations,
  n = n,
  C = C, 
  m = m
)

plot_standardized_MLE(MLEs = res_exp$uncensored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$uncensored_N,
  is_censored = FALSE
)

plot_standardized_MLE(MLEs = res_exp$censored_MLEs, 
  lambda = lambda, 
  n_s = res_exp$censored_N,
  C = C,
  is_censored = TRUE
)

hist(res_exp$censored_N, breaks = 30)
# gets cut off at m, otherwise normally distributed around WHAT VALUE?



# Observation -----------------------------------------
"""
The Fisher Info I(lambda) = D / lambda^2, whereas D is different for 
different censoring schemes. If we choose D correctly then the standardized 
MLEs are N(0, 1) distributed (note that hat lambda ~ N(lambda, D / lambda^2))

Uncensored: D = n (length of data vector)
Type I: D = N (Number of observations before C)
Type II: D = m (predetermined parameter/ no randomness)
Hybrid Type I: D = min(m, N) where N is the number of observations before C
Hybrid Type II: D = max(m, N)

So D is always the number of actually observed uncensored observations.
If we devide by n instead of the correct D, then we will still get a normal distn
but it is 
sqrt(n) hat lambda - lambda / lambda  ~  N(0, n / D)
where n / D > 1, so the variance of the standardized MLEs is greater.

Problem: We need to keep track of the number of uncensored observations. 
That is easy in the case of type II censoring (always D = m), but for other 
censoring schemes, D is random. 
"""