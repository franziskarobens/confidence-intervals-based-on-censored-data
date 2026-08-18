# Replicate the Figures 1a and 1b in Sundberg
# For this get the multipliers c_l and c_u from all methods (C1-C7) proposed in Sundberg
# Plot the multipliers as a function of the uncensored observations N
# To replicate the exact Figures 1a and 1b use the ball_bearing_data and alpha = 0.05

library(ggplot2)
source("Sundberg/CI_Sundberg.R")

plot_CI_limits <- function(t, alpha = 0.05) {
  N_range <- 1:length(t)
  
  methods <- list(
    C1_Sundberg    = C1_Sundberg,
    C2_Sundberg    = C2_Sundberg,
    C3_Sundberg    = C3_Sundberg,
    C4_Sundberg    = C4_Sundberg,
    C5_Sundberg    = C5_Sundberg,
    C6_Sundberg    = C6_Sundberg,
    C7_Sundberg    = C7_Sundberg
  )

  #
  results <- list()

  for (N in N_range) {
    t_c <- sort(t)[N]
    MLE <- MLE_Sundberg(t, t_c)
    
    for (method_name in names(methods)) {
      CI <- tryCatch(
        methods[[method_name]](t, t_c, alpha = alpha),
        error = function(e) {
          cat("ERROR in", method_name, "at N =", N, ":", conditionMessage(e), "\n")
          c(NA_real_, NA_real_)
        }
      )
      results[[length(results) + 1]] <- data.frame(
        N = N,
        method = method_name,
        CI_l = CI[1], # lower limit of the CI
        CI_u = CI[2], # upper limit of the CI
        c_l = CI[1] / MLE, # lower multiplier
        c_u = CI[2] / MLE # upper multiplier
      )
    }
  }

  # transform list of dataframes in one big dataframe
  results_df <- do.call(rbind, results)

  # plot lower c_l
  ggplot(results_df, aes(x = N, y = c_l, color = method)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 0.8, na.rm = TRUE) +
    labs(
      x = "N (number of uncensored observations)",
      y = expression(c[l]),
      title = "Lower CI limit vs. N",
      color = "Method"
    ) +
    theme_minimal()

  # plot upper c_u
  ggplot(results_df, aes(x = N, y = c_u, color = method)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 0.8, na.rm = TRUE) +
    labs(
      x = "N (number of uncensored observations)",
      y = expression(c[l]),
      title = "Upper CI limit vs. N",
      color = "Method"
    ) +
    theme_minimal()

}

# Get Figures 1a and 1b from Sundberg
# plot_CI_limits(t = ball_bearing_data, alpha = 0.05)

