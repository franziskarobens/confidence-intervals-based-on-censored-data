# Parametric Bootstrap-t Procedure

# Simulate the total time on test T = sum(Y_i) of one type II 
# exponentially distributed censored sample
sim_T <- function(n, r, theta) {
  x <- sort(rexp(n, rate = theta))
  sum(x[1:r]) + (n - r) * x[r]          # = sum(pmin(x, x[r]))
}

# PBT Procedure for Exp distributed RV
PBT_CI <- function(T_obs, n, r, alpha = 0.1, B = 2000) {

  theta_hat <- r / T_obs                 # MLE
  se_hat    <- theta_hat / sqrt(r)       # sqrt(1 / I(theta_hat)), I = r / theta^2

  # Bootstrap-t: resample from the fitted model Exp(theta_hat)
  z_star <- replicate(B, {
    T_b     <- sim_T(n, r, theta_hat)
    theta_b <- r / T_b                   # bootstrap MLE
    se_b    <- theta_b / sqrt(r)         # bootstrap version of se
    (theta_b - theta_hat) / se_b
  })

  q <- quantile(z_star, probs = c(1 - alpha / 2, alpha / 2), names = FALSE)

  c(lower = theta_hat - q[1] * se_hat,
    upper = theta_hat - q[2] * se_hat)
}

# Example
n <- 30
r <- 10
x <- sort(rexp(n, rate = 50))
T_obs <- sum(x[1:r]) + (n - r) * x[r]

PBT_CI(T_obs, n = n, r = r, alpha = 0.1, B = 2000)

