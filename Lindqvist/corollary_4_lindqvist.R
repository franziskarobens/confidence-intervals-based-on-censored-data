# Lindqvist Section 5

n <- 2
U <- rexp(n = n, rate = 1)

c <- 1

theta <- 2


compute_tau <- function(theta, U, c) {

  numerator <- sum(U / theta < c)
  denominator <- sum(pmin(U / theta, c))

  return(numerator / denominator)
}

plot_figure_1 <- function(U, c, theta_range = c(0.01, 10), n_points = 500) {
  thetas <- seq(theta_range[1], theta_range[2], length.out = n_points)
  tau_values <- sapply(thetas, compute_tau, U = U, c = c)

  plot(thetas, tau_values, type = "l", col = "blue", lwd = 2,
       xlab = expression(theta), 
       ylab = expression(tau(bold(U), theta))
      )
  abline(v = U / c, col = "gray50", lty = 3)
}

plot_figure_1(U = U, c = c, theta_range = c(0.01, 1.3))

solve_theta_for_u <- function(U, t, c) {
  # use (A1) in the Appendix A from Lindqvist

  sorted_U <- sort(U)
  n <- length(U)

  green_bound <- 1 / (n * c)
  yellow_bounds <- numeric(n - 1) 
  orange_bounds <- numeric(n) 

  for (i in 1:n) {
    numerator <- i * sorted_U[i]
    denominator <- sum(sorted_U[1:i]) + (n - i) * sorted_U[i]
    denominator <- denominator * c

    orange_bounds[i] <- numerator / denominator
  }

  for (i in 1:(n - 1)) {
    numerator <- i * sorted_U[i + 1]
    denominator <- sum(sorted_U[1 : (i+1)]) + (n - i - 1) * sorted_U[i+1]
    denominator <- denominator * c

    yellow_bounds[i] <- numerator / denominator
  }

  pink_bounds <- c(NA, yellow_bounds)

  # case A
  if (t >= 0 & t <= green_bound) {
    cat("Case A\n")
    return(sorted_U[1] / c)
  }

  # Case B
  for (i in 2:n) {
    if (t > pink_bounds[i] & t <= orange_bounds[i]) {
      cat("Case B\n")
      return(sorted_U[i] / c)
    }
  }

  # Case C
  for (i in 1:(n-1)) {
    if (t > orange_bounds[i] & t <= yellow_bounds[i]) {
      cat("Case C\n")
      return(t * sum(sorted_U[1:i]) / (i - (n - i) * c * t))
    }
  }  

  # Case D
  if (t > orange_bounds[n]) {
    cat("Case D\n")
    return(t * sum(sorted_U) / n)
  }

  stop("No result for solve_theta_for_u()")

}


# ----- Inference based on asymptotic distribution ------------





MLE_asymptotic_inference <- function(theta, U, c) {

  X <- U / theta
  Y <- pmin(X, c)

  R <- sum(Y == X)
  S <- sum(Y)

  MLE <- R / S
  return(MLE) # should be theta, works only if n big
}

# -------- Use Corollary 4 from Lindqvist --------------

CI_with_Corollary_4 <- function(theta, U, c, k) {

  n <- length(U)
  T <- compute_tau(theta, U, c)

  if (T == 0) {
    stop("We condition on t > 0 to use Corollary 4.")
  }

  # draw m  samples U
  m <- 100      # large number

  bar_theta <- numeric(m)

  for (i in 1:m) {
    U <- rexp(n = n, rate = 1)
    # solve for theta 
    bar_theta[i] <- solve_theta_for_u(U, t = T, c)
  }

  sorted_theta <- sort(bar_theta)

  cat("Lower bound for theta with alpha = ", round(1 - k/(m+1), 3), "   :   ", sorted_theta[k], "\n")
  cat("Upper bound for theta with alpha = ", round(1 - k/(m+1), 3), "   :   ", sorted_theta[m-k+1], "\n")
  cat("CI for theta with alpha = ", round(1 - 2 * k/(m+1), 3), "   :   ",
      "[", sorted_theta[k], ",", sorted_theta[m - k + 1], "]", "\n")
}

theta <- 2
n <- 10
U <- rexp(n = n, rate = 1)
c <- 1

CI_with_Corollary_4(theta, U, c, k = 5)

