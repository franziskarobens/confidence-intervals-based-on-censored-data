# Replicate Figures 2-5 in Sundberg
# Use calculation Methods presented in Section 3 of Sundberg

library(ggplot2)
source("Sundberg/truncated_exp_distn.R")

theta <- 1
C <- theta / 2
n_range <- 4:50
alpha <- 0.1



n <- 4


# number of uncensored observations
N <- rbinom(n = 1, size = n, prob = 1 - exp(- C / theta))


censored_sample <- rep(C, n - N)
uncensored_sample <- rtexp(n = N, rate = 1, b = C)
sample <- c(censored_sample, uncensored_sample)

# coverage probability
theta_hat <- sum(sample) / N


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


# ---- 1. Characteristic function of the truncated exponential ----

char_trunc_exp <- function(t, C, theta) {
  # t can be a numeric vector; returns complex vector
  numerator   <- 1 - exp(-C * (1 - 1i * theta * t) / theta)
  denominator <- 1 - exp(-C / theta)
  (1 / (1 - 1i * theta * t)) * numerator / denominator
}


# ---- 2. CDF of the sum of N iid truncated-exponential r.v.s ----
#         via the Fourier inversion formula from the paper

F_sum_given_N <- function(z, N, C, theta, Tmax = 300, subdivisions = 3000L) {

  if (N == 0) {
    # Sum of zero terms is degenerate at 0
    return(as.numeric(z >= 0))
  }

  # Since phi(-t) = Conj(phi(t)) and sin(zt/2)/t is even in t,
  # the integrand satisfies integrand(-t) = Conj(integrand(t)),
  # so the integral over [-T, T] equals 2 * Int_0^T Re(integrand(t)) dt,
  # and F(z|N) = (2/pi) * Int_0^T Re(integrand(t)) dt.

  integrand_re <- function(t) {
    phi_pow <- char_trunc_exp(t, C, theta)^N
    val <- (sin(z * t / 2) / t) * exp(-1i * t * z / 2) * phi_pow
    Re(val)
  }

  # Remove the removable singularity at t = 0
  # (lim_{t->0} sin(zt/2)/t = z/2, and phi(0)^N = 1)
  integrand_re_safe <- Vectorize(function(t) {
    if (abs(t) < 1e-10) return(z / 2)
    integrand_re(t)
  })

  I <- stats::integrate(
    integrand_re_safe,
    lower = 0, upper = Tmax,
    subdivisions = subdivisions,
    rel.tol = 1e-8,
    stop.on.error = FALSE
  )$value

  Fz <- (2 / pi) * I
  # Numerical guard: clip tiny overshoots from oscillatory integration
  # min(max(Fz, 0), 1)
  return(Fz)
}


# ---- 3. Exact coverage probability P_theta(theta_hat < t) ----

coverage_prob <- function(t, n, C, theta, Tmax = 300) {

  p <- 1 - exp(-C / theta)   # P(a single observation is uncensored)

  total <- 0
  for (N in 0:n) {
    bin_prob <- dbinom(N, size = n, prob = p)
    # if (bin_prob < 1e-14) next   # skip negligible terms (speed)

    z  <- N * t - (n - N) * C
    Fz <- F_sum_given_N(z, N, C, theta, Tmax = Tmax)

    total <- total + bin_prob * Fz
  }

  total
}


# ---- 4. non coverage probability for all Ci ----

# Settings

n_range <- 4:50
theta <- 1
alpha <- 0.01
C <- theta / 2
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

pb <- txtProgressBar(
  min = 0,
  max = length(n_range) * 7,
  style = 3
)

progress <- 0

# Main calculation

for (n in n_range) {

  # Generate one sample for this n.
  # The SAME sample is used for all 7 methods.
  y <- rexp(
    n = n,
    rate = 1 / theta
  )

  # If all observations are censored,
  # coverage_prob() cannot be used.
  if (all(y >= C)) {
    progress <- progress + n_methods
    setTxtProgressBar(pb, progress)
    next
  }

  # Calculate C1_Sundberg() ... C7_Sundberg()
  for (i in 1:7) {

    # Get the appropriate function:
    # C1_Sundberg, C2_Sundberg, ..., C7_Sundberg
    CI_fun <- get(
      paste0("C", i, "_Sundberg")
    )

    # Calculate confidence interval
    CI <- CI_fun(
      t = y,
      t_c = C,
      alpha = alpha
    )

    # Lower endpoint
    p_lower <- coverage_prob(
      t = CI[1],
      n = n,
      C = C,
      theta = theta,
      Tmax = Tmax
    )

    # Upper endpoint
    p_upper <- coverage_prob(
      t = CI[2],
      n = n,
      C = C,
      theta = theta,
      Tmax = Tmax
    )

    # Noncoverage probability
    non_coverage_prob[
      as.character(n),
      i
    ] <- 1 + p_lower - p_upper

    cat("n = ", n, ", i = ", i, " result = ", 1 + p_lower - p_upper, "\n")

    # Update progress bar
    progress <- progress + 1
    #setTxtProgressBar(pb, progress)
  }
}

close(pb)

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
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right"
  ) +
  geom_hline(yintercept = alpha, color = "black", linetype = "dashed", linewidth = 0.5)

# TODO: die Simulation zeigt noch nicht die gewünschte Non Coverage Probability. 
# entweder sind die CI falsch oder die simulation an sich (wahrscheinlicher)