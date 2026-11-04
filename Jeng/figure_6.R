n <- 50
p_f <- 0.1
t_c <- - log(1 - p_f)
n_sim <- 100

library(ggplot2)
library(survival)
library(scales)

set.seed(2026)

# Setup
mu <- 0
sigma <- 1
n  <- 50
p_f <- 0.1
t_c <- -log(1 - p_f)          # censoring time, ~0.105
nsim <- 100

# Time grid
# t_grid <- exp(seq(log(0.005), log(2), length.out = 200))
t_grid <- seq(0.005, 2, length.out = 200)


# Simulate and fit
fits <- list()
for (i in 1:nsim) {
  x <- rweibull(n, shape = 1 / sigma, scale = exp(mu))
  status <- as.numeric(x <= t_c)
  time   <- pmin(x, t_c)

  if (sum(status) < 2) next            # ML estimate of (mu, sigma) needs at least 2 failures

  fit <- tryCatch(
    survreg(Surv(time, status) ~ 1, dist = "weibull"),
    error = function(e) NULL
  )
  if (is.null(fit)) next

  mu_hat    <- unname(coef(fit))
  sigma_hat <- fit$scale
  fits[[length(fits) + 1]] <- data.frame(
    sim = length(fits) + 1,
    t   = t_grid,
    F   = pweibull(t_grid, shape = 1 / sigma_hat, scale = exp(mu_hat))
  )
}
df <- do.call(rbind, fits)

# True cdf
df_true <- data.frame(t = t_grid,
                      F = pweibull(t_grid, shape = 1 / sigma, scale = exp(mu)))

# Weibull probability scale: y = log(-log(1 - p))
sev_trans <- trans_new(
  "sev",
  transform = function(p) log(-log(1 - p)),
  inverse   = function(z) 1 - exp(-exp(z))
)

ggplot() +
  geom_line(data = df, aes(t, F, group = sim),
            colour = "grey40", alpha = 0.5, linewidth = 0.4) +
  geom_line(data = df_true, aes(t, F),
            colour = "red", linewidth = 1) +
  geom_vline(xintercept = t_c, linetype = "dashed") +
  scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1)) +
  scale_y_continuous(trans = sev_trans,
                     breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.9, 0.99),
                     limits = c(0.001, 0.99)) +
  labs(x = "t (log scale)", y = "Proportion failing",
       title = "Weibull life test simulation: 100 ML estimates of F(t)",
       subtitle = paste0("n = 50, type I censoring at t_c = ", round(t_c, 3),
                         " (p_f = 0.1, E(r) = 5); true F in red"))
