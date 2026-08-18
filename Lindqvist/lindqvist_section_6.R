# --------------------------------------------------------------------------
# Fiducial inference simulation for n = 2 (Section 6 example)
# 
# Setting: c = 1, one failure at Y1 = 0.5, one censoring at Y2 = 1
# => R = 1, S = 3/2, T = t = 2/3
#
# For a simulated pair u = (u1, u2) with order statistics u(1) <= u(2):
#
#   bar theta(u, t) = 2*u(1)   if u(2) >= 2*u(1)
#   bar theta(u, t) = u(2)     if u(2) <  2*u(1)
#
# Fiducial distribution (5):
#   theta(u, t) = 2*u(1)   conditional on  u(2) >= 2*u(1)
# --------------------------------------------------------------------------

set.seed(123)

c_val <- 1
t_val <- 2 / 3
n_sim <- 100000

# --- Simulate 100,000 pairs of iid standard exponential variables ---------
U1 <- rexp(n_sim, rate = 1)
U2 <- rexp(n_sim, rate = 1)

u_min <- pmin(U1, U2)   # u(1)
u_max <- pmax(U1, U2)   # u(2)

# --- theta(u, t): unconditional piecewise solution -------------------------
# theta = 2*u(1) if u(2) >= 2*u(1), else theta = u(2)
theta_full <- ifelse(u_max >= 2 * u_min, 2 * u_min, u_max)

# --- Fiducial version (5): only keep draws where u(2) >= 2*u(1) -----------
fiducial_condition <- u_max >= 2 * u_min
theta_fiducial <- 2 * u_min[fiducial_condition]

cat("Fraction of draws satisfying u(2) >= 2*u(1):",
    round(mean(fiducial_condition), 4), "\n\n")

# --- Confidence intervals: delete 2.5% at each end of the histogram -------
CI_full <- quantile(theta_full, probs = c(0.025, 0.975))
CI_fiducial <- quantile(theta_fiducial, probs = c(0.025, 0.975))

cat("95% CI (theta(u, 2/3), unconditional):  [",
    round(CI_full[1], 3), ",", round(CI_full[2], 3), "]\n")
cat("95% CI (fiducial version, conditional): [",
    round(CI_fiducial[1], 3), ",", round(CI_fiducial[2], 3), "]\n\n") 
# -> same values as in paper

# --------------------------------------------------------------------------
# Plot: two histograms side by side, matching "Figure 2" in the text
# --------------------------------------------------------------------------
par(mfrow = c(1, 2))

hist(theta_full, breaks = 100, freq = FALSE, xlim = c(0, 5),
     main = expression(theta(u, t == 2/3)),
     xlab = expression(theta), ylab = "Density",
     col = "lightblue", border = "white")
abline(v = CI_full, col = "red", lty = 2, lwd = 2)

# Given: F_theta(x) = 2*exp(-2*x) - exp(-3*x)   (survival function)
# Density = -d/dx F(x) = 4*exp(-2*x) - 3*exp(-3*x)

x_grid <- seq(0, 5, length.out = 500)
density_theta <- 4 * exp(-2 * x_grid) - 3 * exp(-3 * x_grid)

lines(x_grid, density_theta, col = "darkred", lwd = 2)

legend("topright", legend = "Theoretical density", col = "darkred", lwd = 2, bty = "n")

hist(theta_fiducial, breaks = 100, freq = FALSE, xlim = c(0, 5),
     main = "Fiducial distribution (5)",
     xlab = expression(theta), ylab = "Frequency",
     col = "lightgreen", border = "white")
abline(v = CI_fiducial, col = "red", lty = 2, lwd = 2)

x_grid <- seq(0, 5, length.out = 500)
density_theta <- 3 / 2 * exp(- 3/2 * x_grid)

lines(x_grid, density_theta, col = "darkred", lwd = 2)

legend("topright", legend = "Theoretical density Exp(3/2)", col = "darkred", lwd = 2, bty = "n")


par(mfrow = c(1, 1))


# ----- Helper function to get theta using eq (4) --------------------------------

theta_from_u <- function(u, t, c) {

  n <- length(u)

  u <- sort(u)

  A <- c(0, cumsum(u))

  theta_best <- 0

  for (k in 0:n) {

    if (k == 0) {
      lower <- 0
    } else {
      lower <- u[k] / c
    }

    if (k == n) {
      upper <- Inf
    } else {
      upper <- u[k + 1] / c
    }

    if (k == 0) {

      theta_best <- max(theta_best, upper)

    } else {

      A_k <- A[k + 1]

      if (k == n) {

        candidate <- t * A_k / n

        if (candidate >= lower) {
          theta_best <- max(theta_best, candidate)
        }


      } else {
        denominator <- k - t * (n - k) * c

        tau_upper <- k * upper /
          (A_k + (n - k) * c * upper)

        if (t >= tau_upper) {

          theta_best <- max(theta_best, upper)

        } else {

          if (denominator > 0) {

            candidate <- t * A_k / denominator

            if (candidate >= lower && candidate <= upper) {
              theta_best <- max(theta_best, candidate)
            }
          }
        }
      }
    }
  }

  return(theta_best)
}


# ----- Check coverage probabilities of bounds and CI -----------------------------------------

n <- 10
c <- 1
theta_true <- 1
n_simulation <- 10000
m <- 1000

R_values <- numeric(n_simulations)
T_values <- numeric(n_simulations)
lower_bounds <- numeric(n_simulations)
upper_bounds <- numeric(n_simulations) 

pb <- txtProgressBar(
  min = 0,
  max = n_simulation,
  style = 3
)

for (sim in 1:n_simulation) {

  # Generate n = 10 exponential lifetimes with true hazard theta = 1
  X <- rexp(n, rate = theta_true)

  # Type-I censoring at c = 1
  Y <- pmin(X, c_val)


  # Number of observed failures
  R <- sum(X < c_val)
  R_values[sim] <- R

  # Total observed time
  S <- sum(Y)


  # MLE
  # T = R/S
  t <- R / S
  T_values[sim] <- t


  # Generate m = 1000 vectors U
  # Each row is one vector: U = (U1,...,U10) where Ui ~ Exp(1)
  U_matrix <- matrix(
    rexp(m * n, rate = 1),
    nrow = m,
    ncol = n
  )


  # Calculate theta(U,t) for each of the m simulated U-vectors
  theta_values <- numeric(m)

  for (j in 1:m) {

    theta_values[j] <- theta_from_u(
      U_matrix[j, ],
      t,
      c_val
    )
  }


  # Fiducial/confidence limits: remove 2.5% from each tail
  bounds <- quantile(
    theta_values,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE
  )


  lower_bounds[sim] <- bounds[1]
  upper_bounds[sim] <- bounds[2]


  # Progress indicator
  setTxtProgressBar(pb, sim)
}

close(pb)


# Coverage calculations
lower_coverage <- mean(lower_bounds <= theta_true)

upper_coverage <- mean(upper_bounds >= theta_true)

two_sided_coverage <- mean(
  lower_bounds <= theta_true &
  upper_bounds >= theta_true
)


# Results

cat("\n")
cat("============================================================\n")
cat("Simulation results\n")
cat("============================================================\n")

cat("n =", n, "\n")
cat("c =", c_val, "\n")
cat("true theta =", theta_true, "\n")
cat("number of datasets =", n_sim, "\n")
cat("m =", m, "\n\n")

cat("Coverage of lower bound:\n")
cat(round(100 * lower_coverage, 2), "%\n\n")

cat("Coverage of upper bound:\n")
cat(round(100 * upper_coverage, 2), "%\n\n")

cat("Coverage of two-sided interval:\n")
cat(round(100 * two_sided_coverage, 2), "%\n\n")

cat("Nominal coverage = 95%\n")



# ------ 