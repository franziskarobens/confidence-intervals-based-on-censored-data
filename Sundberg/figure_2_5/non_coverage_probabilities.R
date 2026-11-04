# Replicate Figures 2-5 in Sundberg
# Use calculation Methods presented in Section 3 of Sundberg

library(ggplot2)
library(reticulate)
py_install("scipy")   # installs into reticulate's active Python env
source("Sundberg/truncated_exp_distn.R")
source_python("Sundberg/calc_integral_F.py")


# ============================================================
# Exact coverage probabilities for Type I censored exponential
# data, following the Fourier-inversion approach (Sundberg-type
# method) described in the paper excerpt.
#
# Model:
#   - theta = mean of the exponential distribution (rate = 1/theta)
#   - n     = sample size
#   - C     = fixed censoring time (Type I censoring)
#   - N     = number of uncensored observations,
#             N ~ Binomial(n, p),  p = 1 - exp(-C/theta)
#   - Given N, the N uncensored values Z_1,...,Z_N are iid draws
#     from Exp(theta) truncated to [0, C]
#   - theta_hat = (1/N) * sum_{i=1}^n Y_i   (Y_i = observed,
#     possibly censored, value)
#
# Goal: compute P_theta(theta_hat < t) exactly via
#   P_theta(theta_hat < t)
#     = E_N[ P_theta( sum_{i=1}^N Z_i < N*t - (n-N)*C | N ) ]
#
# where the inner probability is obtained by inverting the
# characteristic function of the truncated exponential via the
# Levy/Gil-Pelaez-type inversion formula given in the paper:
#
#   F(z|N) = lim_{T->inf} (1/pi) * Int_{-T}^{T}
#                [sin(z t/2)/t] * exp(-i t z/2) * phi(t)^N  dt
#
# with phi(t) the characteristic function of the truncated
# exponential (truncated at C, mean theta):
#
#   phi(t) = 1/(1 - i*theta*t) *
#            (1 - exp(-C*(1 - i*theta*t)/theta)) /
#            (1 - exp(-C/theta))
# ============================================================


# ---- Exact coverage probability P_theta(theta_hat < t) ----

coverage_prob <- function(t, n, C, theta, Tmax = 300) {

  p <- 1 - exp(- C / theta)   # P(a single observation is uncensored)

  total <- 0
  for (N in 0:n) {
    bin_prob <- dbinom(N, size = n, prob = p)
    if (bin_prob < 1e-14) next   # skip negligible terms (speed)

    z  <- N * t - (n - N) * C
    Fz <- fourier_inversion_integral(z = z, N = N, C = C, theta = theta, T = Tmax)

    total <- total + bin_prob * Fz
  }

  total
}


# ---- Non coverage probability for all CI (Simulation) ----
# Use the method proposed in Section 3 of Sundberg

# Settings
n_range <- 4:50
theta <- 1
alpha <- 0.1
C <- theta
Tmax <- 300

# Storage for noncoverage probabilities
non_coverage_prob <- matrix(
  NA_real_,
  nrow = length(n_range),
  ncol = 7,
  dimnames = list(
    n_range,
    paste0("C", 1:7)
  )
)



# Progress bar
pb <- txtProgressBar(min = 0, max = length(n_range) * 7, style = 3)
progress <- 0

# Main calculation
for (n in n_range) {

  # Calculate C1_Sundberg() ... C7_Sundberg()
  for (i in 1:7) {

    if(i == 5) next
    
    CI_fun <- get(paste0("C", i, "_Sundberg"))

    # Calculate confidence interval
    CI <- CI_fun(N = N, t_c = C, alpha = alpha)
    c_lower <- CI$c_lower
    c_upper <- CI$c_upper

    # Lower endpoint
    p_lower <- 0
    for (j in 1:n) {
      c_lower <- 
      coverage_prob(t = theta / c_lower, n = n, C = C, theta = theta, Tmax = Tmax)
    }

    p_lower <- coverage_prob(t = theta / c_lower, n = n, C = C, theta = theta, Tmax = Tmax)

    # Upper endpoint
    p_upper <- coverage_prob(t = theta / c_upper, n = n, C = C, theta = theta, Tmax = Tmax)

    # Noncoverage probability
    non_coverage_prob[as.character(n), i] <- 1 - p_lower + p_upper

    # Update progress bar
    progress <- progress + 1
    setTxtProgressBar(pb, progress)
  }
}

close(pb)


calc_CI_bounds <- function(k, alpha) {
  CI_bounds <- data.frame(
    c_l = numeric(7),
    c_u = numeric(7),
    row.names = as.character(1:7)
  )

  for (i in 1:7) {
    if (i == 5) next
    CI_fun <- get(paste0("C", i, "_Sundberg"))
    CI <- CI_fun(N = k, alpha = alpha)
    CI_bounds[as.character(i), "c_l"] <- CI$c_l
    CI_bounds[as.character(i), "c_u"] <- CI$c_u
  }
  return (CI_bounds)
}

# Plot the CI width vs N
k_range <- 0:10

width_df <- do.call(rbind, lapply(k_range, function(k) {
  CI_bounds <- calc_CI_bounds(k = k, alpha = 0.1)
  data.frame(
    k = k,
    method = paste0("C", 1:7),
    width = CI_bounds$c_u - CI_bounds$c_l
  )
}))

ggplot(width_df, aes(x = k, y = width, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.8) +
  labs(
    x = "N (number of uncensored observations)",
    y = expression(c[u] - c[l]),
    title = "CI width vs. k",
    color = "Method"
  ) +
  theme_minimal()
# ENDE

non_coverage_prob <- as.data.frame(
  matrix(NA_real_, nrow = length(n_range), ncol = 7,
         dimnames = list(as.character(n_range), paste0("C", 1:7)))
)

z <- qnorm(p = 1 - alpha / 2)
# Only for C1
for (n in n_range) {

  p_lower <- numeric(7) # P(hat theta < theta / c_l)
  p_upper <- numeric(7) # P(hat theta < theta / c_u)

  for (k in 2:n) { # hier muss eigentlich k von 0 starten!!
    p <- dbinom(k, size = n, prob = 1 - exp(- C / theta))

    CI_bounds <- calc_CI_bounds(k, alpha)
    arguments <- k * theta / CI_bounds - (n - k) * C
    
    Fz <- as.data.frame(
      lapply(arguments, function(col) {
        sapply(col, function(z) {
          fourier_inversion_integral(T = 100, z = z, N = k, theta = theta, C = C)
        })
      })
    ) 
    # 1) Problem here: Fz is negative in many cases, causing that Fz = 0
    # 2) Another Problem: k = 0, then CI_bounds doesn't work and p_... becomes NaN
    # 3) Also for k = 0, 1, bounds of C2 dont fulfil c_l < c_u
    
    p_lower <- p_lower + p * Fz[, 1]
    p_upper <- p_upper + p * Fz[, 2]
    cat("n = ", n,
        ", k = ", k, "\n",
        "Fz = ", Fz[, 1], "\n\n")
  }

  non_coverage_prob[as.character(n), ] <- 1 - p_lower + p_upper
}

# Convert results to data frame for ggplot

results <- as.data.frame(non_coverage_prob)

results$n <- n_range

results_long <- reshape(
  results,
  varying = paste0("C", 1:7),
  v.names = "noncoverage",
  timevar = "method",
  times = paste0("C", 1:7),
  direction = "long"
)

# Clean row names
rownames(results_long) <- NULL

# Plot

title_text <- paste(
  "Noncoverage probabilities for Sundberg CIs with",
  "theta =", theta,
  ", C =", C,
  "and alpha =", alpha
)

ggplot(
  results_long,
  aes(
    x = n,
    y = noncoverage,
    color = method
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(
    title = title_text,
    x = "Sample size n",
    y = "Noncoverage probability",
    color = "Method"
  ) +
  ylim(0.06, 0.14) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right"
  ) +
  geom_hline(yintercept = alpha, color = "black", linetype = "dashed", linewidth = 0.5)

# TODO: die Simulation zeigt noch nicht die gewünschte Non Coverage Probability. 
# entweder sind die CI falsch oder die simulation an sich (wahrscheinlicher)