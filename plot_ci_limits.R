# this file is for plotting the limits of the CI based on how many observations are censored
# see Figures 1a and 1b in Sundberg


library(ggplot2)

plot_ci_limits <- function(t, alpha) {
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

  results <- list()

  for (N in N_range) {
    t_c <- sort(t)[N]
    mle <- mle_Sundberg(t, t_c)
    
    for (method_name in names(methods)) {
      ci <- tryCatch(
        methods[[method_name]](t, t_c, alpha = alpha),
        error = function(e) {
          cat("ERROR in", method_name, "at N =", N, ":", conditionMessage(e), "\n")
          c(NA_real_, NA_real_)
        }
      )
      results[[length(results) + 1]] <- data.frame(
        N = N,
        method = method_name,
        ci_l = ci[1],
        ci_u = ci[2], 
        c_l = ci[1] / mle,
        c_u = ci[2] / mle 
      )
    }
  }

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


plot_ci_limits(t = ball_bearing_data, alpha = 0.05)

