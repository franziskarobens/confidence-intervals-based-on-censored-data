# Discover the problem that for C2 in Sundberg, for small n 
# the interval bounds don't fulfil c_l < c_u
# Plot sample size n vs CI width

library(ggplot2)

n <- 0:10

alpha <- 0.1

# Multipliers for C2 from Sundberg
multipliers <- C2_Sundberg(N = n, alpha = alpha)

c_upper <- multipliers$c_upper
c_lower <- multipliers$c_lower

widths <- c_upper - c_lower

plot_df <- data.frame(n = n, width = widths)

ggplot(plot_df, aes(x = n, y = width)) +
  geom_point(size = 2) +
  labs(
    x = "n",
    y = expression(c[u] - c[l]),
  ) +
  scale_x_continuous(breaks = n) +
  theme_minimal()

# For all n < z, it holds c_l > c_u.
qnorm(1 - alpha / 2)^2

