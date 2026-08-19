# Experiment: Compare exact and asymptotic confidence intervals
# The non coverage probability of the asymptotic CI should converge to alpha
# but for the exact CI, it should be alpha for all n (even for small)
# In this example we use CIs for the mean of a normal distribution

library(ggplot2)

mu <- 10
sigma = 100
alpha <- 0.1
z <- qnorm(1 - alpha / 2)

number_simulations <- 1000
n_max <- 100

non_coverage_prob_exact <- numeric(n_max-1)
non_coverage_prob_asymptotic <- numeric(n_max-1)

for (n in 2:n_max) {

  not_in_exact_CI <- 0
  not_in_asymptotic_CI <- 0

  for (i in 1:number_simulations) {
    X <- rnorm(n = n, mean = mu, sd = sigma)
    X_bar <- mean(X)
    S <- sd(X)

    exact_CI <- c(
      X_bar - z * sigma / sqrt(n), 
      X_bar + z * sigma / sqrt(n)
    )


    asymptotic_CI <- c(
      X_bar - z * S / sqrt(n), 
      X_bar + z * S / sqrt(n)
    )

    if (mu < exact_CI[1] | mu > exact_CI[2]) {
      not_in_exact_CI <- not_in_exact_CI + 1
    }

    if (mu < asymptotic_CI[1] | mu > asymptotic_CI[2]) {
      not_in_asymptotic_CI <- not_in_asymptotic_CI + 1
    }
  }
  
  non_coverage_prob_exact[n-1] <- not_in_exact_CI / number_simulations
  non_coverage_prob_asymptotic[n-1] <- not_in_asymptotic_CI / number_simulations
}


results_df <- data.frame(
  n = rep(2:n_max, times = 2),
  non_coverage = c(non_coverage_prob_exact, non_coverage_prob_asymptotic),
  type = rep(c("Exact CI", "Asymptotic CI"), each = n_max - 1)
)

ggplot(results_df, aes(x = n, y = non_coverage, color = type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("Exact CI" = "red", "Asymptotic CI" = "blue")) +
  labs(
    x = "n (sample size)",
    y = "P(non-coverage)",
    title = "Non-coverage probability: exact vs. asymptotic CI",
    color = "CI type"
  ) +
  theme_minimal()
