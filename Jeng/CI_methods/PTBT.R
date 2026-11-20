# Parametric Transformed Bootstrap Procedure

PTBT_CI <- function(T_obs, n, r, g, alpha = 0.1, B = 2000) {

  theta_hat <- r / T_obs
  se_g_hat  <- abs(g$g_der(theta_hat)) * theta_hat / sqrt(r)

  # Bootstrap: resample from the fitted model Exp(theta_hat)
  z_star <- replicate(B, {
    T_b     <- sim_T(n, r, theta_hat)
    theta_b <- r / T_b                              # bootstrap MLE
    se_g_b  <- abs(g$g_der(theta_b)) * theta_b / sqrt(r) # bootstrap version of se_g
    (g$g(theta_b) - g$g(theta_hat)) / se_g_b
  })

  q <- quantile(z_star, probs = c(1 - alpha / 2, alpha / 2), names = FALSE)

  # Endpoints on the g-scale: g(theta_hat) - z_(1-a/2) * se_g and g(theta_hat) - z_(a/2) * se_g
  end_g <- c(g$g(theta_hat) - q[1] * se_g_hat,
             g$g(theta_hat) - q[2] * se_g_hat)

  # Back-transform; g monotone increasing keeps the order, decreasing reverses it
  end <- g$g_inv(end_g)
  end <- sort(end)

  c(lower = end[1], upper = end[2])
}


# Example
n <- 30
r <- 10
x <- sort(rexp(n, rate = 0.5))
T_obs <- sum(x[1:r]) + (n - r) * x[r]

# TODO: PBT und PTBT mit identity stimmen nicht überein!
PBT_CI(T_obs, n, r, B = 2000)
PTBT_CI(T_obs, n, r, g = g_identity(), alpha = 0.1, B = 2000)

PTBT_CI(T_obs, n, r, g = g_exp(), alpha = 0.1, B = 2000)
PTBT_CI(T_obs, n, r, g = g_log(), alpha = 0.1, B = 2000)
PTBT_CI(T_obs, n, r, g = g_one_third(), alpha = 0.1, B = 2000)
