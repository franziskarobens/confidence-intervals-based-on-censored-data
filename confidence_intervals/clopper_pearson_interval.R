
# --- Example: Clopper-Pearson Interval (Section 2.1 in Lindqvist) --------------

n <- 100
alpha <- 0.05

n_simulations <- 10000000
theta_in_upper_CI <- logical(n_simulations)
theta_in_lower_CI <- logical(n_simulations)

theta <- 0.5

a <- numeric(n+1)
b <- numeric(n+1)

a[n+1] <- 1       # a(n) = 1
b[1] <- 0         # b(0) = 0


for (t in 0:(n - 1)) {
  a[t + 1] <- qbeta(1 - alpha, t + 1, n - t)
}

for (t in 1:n) {
  b[t + 1] <- qbeta(alpha, t, n + 1 - t)
}

T <- rbinom(n = n_simulations, size = n, prob = theta)

# upper CI:
# c(0, a[T+1])

# lower CI:
# c(b[T+1], 1)

theta_in_upper_CI <- theta <= a[T+1]
theta_in_lower_CI <- theta >= b[T+1]



# should be 1 - alpha
sum(theta_in_upper_CI) / length(theta_in_upper_CI)
sum(theta_in_lower_CI) / length(theta_in_lower_CI)

