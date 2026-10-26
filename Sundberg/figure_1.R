# Replicate Figure 1a and 1b in Sundberg

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggh4x)

alpha   <- 0.025
N_range <- 1:25
CI_nrs  <- 1:7

res <- matrix(NA_real_, ncol = 4, nrow = length(N_range) * length(CI_nrs))
colnames(res) <- c("N", "CI_nr", "c_lower", "c_upper")

row <- 1
for (N in N_range) {
  for (nr in CI_nrs) {
    lim <- tryCatch(
      Sundberg_limits(method_nr = nr, r = N, alpha = alpha),
      error = function(e) {
        warning(sprintf("C%d, N = %d: %s", nr, N, conditionMessage(e)))
        list(c_lower = NA_real_, c_upper = NA_real_)
      }
    )
    res[row, ] <- c(N, nr, lim$c_lower, lim$c_upper)
    row <- row + 1
  }
}

res <- as.data.frame(res)
head(res, 14)

# Plot

plot_df <- res |>
  mutate(CI_nr = factor(CI_nr, labels = paste0("C", 1:7))) |>
  pivot_longer(c(c_lower, c_upper),
               names_to = "c", values_to = "value") |>
  mutate(c = factor(c,
                        levels = c("c_lower", "c_upper"),
                        labels = c("c_lower", "c_upper")))


ggplot(plot_df, aes(x = N, y = value, colour = CI_nr, linetype = CI_nr)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ c, nrow = 1, scales = "free_y") +
  facetted_pos_scales(
    y = list(
      c == "c_lower" ~ scale_y_continuous(limits = c(0, 0.8), oob = scales::oob_keep),
      c == "c_upper" ~ scale_y_continuous(limits = c(1, 8), oob = scales::oob_keep)
    )
  ) +
  labs(x = "N", y = "c", colour = "Method", linetype = "Method") +
  theme_bw() +
  theme(legend.position = "bottom")

# hier sind die resultate von C1 und C2 im Vergleich zum Paper unterschiedlich
# weil wir eine andere parametrisierung gewählt haben (rate = theta statt 
# rate = 1 / theta)
