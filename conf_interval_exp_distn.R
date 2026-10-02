# Coverage simulation for confidence intervals for the rate lambda
# of an exponential distribution.
# Pivot: 2 * lambda * sum(X) ~ chi-squared(2n)

set.seed(2026)

n_sim  <- 1e5
n      <- 100
lambda <- 2
alpha  <- 0.05
df     <- 2 * n

# Simulate the sample sums (one per replication) ----------------------------
sum_x <- replicate(n_sim, sum(rexp(n, rate = lambda)))
denom <- 2 * sum_x

# One-sided CI: (0, upper] ---------------------------------------------------
upper_one_sided <- qchisq(1 - alpha, df) / denom
coverage_one_sided <- mean(lambda <= upper_one_sided)

# Two-sided CI: [lower, upper] -----------------------------------------------
lower <- qchisq(alpha / 2, df) / denom
upper <- qchisq(1 - alpha / 2, df) / denom

two_sided <- c(
  left   = mean(lambda < lower),                     # lambda left of the CI
  inside = mean(lambda >= lower & lambda <= upper),  # should be 1 - alpha
  right  = mean(lambda > upper)                      # lambda right of the CI
)

# Results ----------------------------------------------------------------------
cat("Target coverage (1 - alpha):", 1 - alpha, "\n\n")
cat("One-sided CI coverage:", coverage_one_sided, "\n\n")
cat("Two-sided CI:\n")
print(two_sided)
