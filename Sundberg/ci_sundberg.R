# Description
# TBD


# returns MLE given by Sundberg in equation (2) in chapter 2.2
# @ params
# t = c(t_1, ..., t_n) with t_i ~ Exp
# t_c censoring time
MLE_Sundberg <- function(t, t_c) {
  N <- sum(t <= t_c) # number of uncensored
  t_obs <- pmin(t, t_c)

  return(sum(t_obs) / N)
}


# returns CI regarding to Sundberg (see C1)
# @params
# t: is a n dimensional vector of observations whereas we assume 
# they are exponentially distributed
# t_c: censoring time 
# alpha: confidence niveau for the CI
C1_Sundberg <- function(N = NULL, t = NULL, t_c = NULL, alpha = 0.1) {

  if (is.null(t) == is.null(N)) {
    stop("Either 't' OR 'N' must be supplied.")
  }
  if (is.null(N)) {
    N <- sum(t < t_c) # number of uncensored
  }

  z <- qnorm(p = 1 - alpha / 2)
  c_lower <- 1 - z / sqrt(N)
  c_upper <- 1 + z / sqrt(N)

  if (is.null(t)) {
    if (! is.null(t_c)) {
      stop("If t is not given, t_c must also be null (don't need t_c for calculating c_l and c_u).")
    }

    return(list(
      c_lower = c_lower,
      c_upper = c_upper
    ))
  }

  if (is.null(t_c)) {
    stop("t_c must be given, if CI should be calculated")
  }
  
  theta_hat <- MLE_Sundberg(t, t_c)
  C1 <- c(
    theta_hat * c_lower, 
    theta_hat * c_upper
  )
  
  return(list(
      CI = unname(C1),
      c_lower = c_lower,
      c_upper = c_upper
    )
  )
    
}


C2_Sundberg <- function(N = NULL, t = NULL, t_c = NULL, alpha = 0.1) {

  if (is.null(t) == is.null(N)) {
    stop("Either 't' OR 'N' must be supplied.")
  }
  if (is.null(N)) {
    N <- sum(t < t_c) # number of uncensored
  }
   
  z <- qnorm(p = 1 - alpha / 2)
  
  c_lower <- 1 / (1 + z / sqrt(N))
  c_upper <- 1 / (1 - z / sqrt(N))

  if (is.null(t)) {
    if (! is.null(t_c)) {
      stop("If t is not given, t_c must also be null (don't need t_c for calculating c_l and c_u).")
    }

    return(list(
      c_lower = c_lower,
      c_upper = c_upper
    ))
  }
  if (is.null(t_c)) {
    stop("t_c must be given, if CI should be calculated")
  }

  theta_hat <- MLE_Sundberg(t, t_c)
 
  C2 <- c(
    theta_hat * c_lower, 
    theta_hat * c_upper
  )
  
  return(list(
      CI = unname(C2),
      c_lower = c_lower,
      c_upper = c_upper
    )
  )
}


C3_Sundberg <- function(N = NULL, t = NULL, t_c = NULL, alpha = 0.1) {

  if (is.null(t) == is.null(N)) {
    stop("Either 't' OR 'N' must be supplied.")
  }
  if (is.null(N)) {
    N <- sum(t < t_c) # number of uncensored
  }

  z <- qnorm(p = 1 - alpha / 2)
  
  c_lower <- exp(- z / sqrt(N))
  c_upper <- exp(z / sqrt(N))

  if (is.null(t)) {
    if (! is.null(t_c)) {
      stop("If t is not given, t_c must also be null (don't need t_c for calculating c_l and c_u).")
    }

    return(list(
      c_lower = c_lower,
      c_upper = c_upper
    ))
  }

  if (is.null(t_c)) {
    stop("t_c must be given, if CI should be calculated")
  }

  theta_hat <- MLE_Sundberg(t, t_c)

  C3 <- c(
    theta_hat * c_lower, 
    theta_hat * c_upper
  )
  
  return(list(
      CI = unname(C3),
      c_lower = c_lower,
      c_upper = c_upper
    )
  )
}


C4_Sundberg <- function(N = NULL, t = NULL, t_c = NULL, alpha = 0.1) {

  if (is.null(t) == is.null(N)) {
    stop("Either 't' OR 'N' must be supplied.")
  }
  if (is.null(N)) {
    N <- sum(t < t_c) # number of uncensored
  }
  
  z <- qnorm(p = 1 - alpha / 2)
  
  c_lower <- 1 / (1 + z / (3*sqrt(N)))^3
  c_upper <- 1 / (1 - z / (3*sqrt(N)))^3



  if (is.null(t)) {
    if (! is.null(t_c)) {
      stop("If t is not given, t_c must also be null (don't need t_c for calculating c_l and c_u).")
    }

    return(list(
      c_lower = c_lower,
      c_upper = c_upper
    ))
  }

  if (is.null(t_c)) {
    stop("t_c must be given, if CI should be calculated")
  }

  # maximum likelihood estimations
  theta_hat <- MLE_Sundberg(t, t_c)

  C4 <- c(
    theta_hat * c_lower, 
    theta_hat * c_upper
  )
  
  return(list(
      CI = unname(C4),
      c_lower = c_lower,
      c_upper = c_upper
    )
  )
}


C5_Sundberg <- function(t, t_c, alpha = 0.1) {
  theta_hat <- MLE_Sundberg(t, t_c)
  N <- sum(t < t_c) # number of uncensored
  z <- qnorm(p = 1 - alpha / 2)
  
  # h(theta) = 2N * (theta_hat/theta - 1 - log(theta_hat/theta)) - z^2
  # b1, b2 solve h(theta) = 0.
  # h is convex in theta with minimum -z^2 at theta = theta_hat,
  # and h -> +Inf as theta -> 0+ or theta -> Inf, so there are
  # exactly two roots: one below theta_hat, one above.
  h <- function(theta) {
    2 * N * (theta_hat / theta - 1 - log(theta_hat / theta)) - z^2
  }
  
  # lower root: search in (epsilon, theta_hat)
  b1 <- uniroot(h, interval = c(theta_hat * 1e-6, theta_hat))$root
  
  # upper root: search in (theta_hat, large multiple of theta_hat)
  # expand upper bound until a sign change is found, in case theta_hat is small
  upper <- theta_hat * 10
  while (h(upper) < 0) {
    upper <- upper * 10
  }
  b2 <- uniroot(h, interval = c(theta_hat, upper))$root
  
  C5 <- c(b1, b2)
  
  return(list(
      CI = unname(C5)
    )
  )
}



C6_Sundberg <- function(N = NULL, t = NULL, t_c = NULL, alpha = 0.1) {

  if (is.null(t) == is.null(N)) {
    stop("Either 't' OR 'N' must be supplied.")
  }
  if (is.null(N)) {
    N <- sum(t < t_c) # number of uncensored
  }

  q1 <- qchisq(p = 1 - alpha / 2, df = 2 * N)
  q2 <- qchisq(p = alpha / 2, df = 2 * N)
  
  c_lower <- 2 * N / q1
  c_upper <- 2 * N / q2

  if (is.null(t)) {
    if (! is.null(t_c)) {
      stop("If t is not given, t_c must also be null (don't need t_c for calculating c_l and c_u).")
    }

    return(list(
      c_lower = c_lower,
      c_upper = c_upper
    ))
  }

  if (is.null(t_c)) {
    stop("t_c must be given, if CI should be calculated")
  }

  theta_hat <- MLE_Sundberg(t, t_c)

  C6 <- c(
    theta_hat * c_lower, 
    theta_hat * c_upper
  )
  
  return(list(
      CI = unname(C6),
      c_lower = c_lower,
      c_upper = c_upper
    )
  )
}

C7_Sundberg <- function(N = NULL, t = NULL, t_c = NULL, alpha = 0.1) {
  
  if (is.null(t) == is.null(N)) {
    stop("Either 't' OR 'N' must be supplied.")
  }
  if (is.null(N)) {
    N <- sum(t < t_c) # number of uncensored
  }

  q1 <- qchisq(p = 1 - alpha / 2, df = 2 * N + 1)
  q2 <- qchisq(p = alpha / 2, df = 2 * N + 1)

  c_lower <- 2 * N / q1
  c_upper <- 2 * N / q2
  
  if (is.null(t)) {
    if (! is.null(t_c)) {
      stop("If t is not given, t_c must also be null (don't need t_c for calculating c_l and c_u).")
    }
    
    return(list(
      c_lower = c_lower,
      c_upper = c_upper
    ))
  }

  if (is.null(t_c)) {
    stop("t_c must be given, if CI should be calculated")
  }

  theta_hat <- MLE_Sundberg(t, t_c)

  C7 <- c(
    theta_hat * c_lower,
    theta_hat * c_upper
  )
  
  return(list(
      CI = unname(C7),
      c_lower = c_lower,
      c_upper = c_upper
    )
  )
}


# Observation: confidence intervals C1 and C2 differ a lot from the remaining intervals
# C3 - C7 are relatively similar

# for edge cases:
# C = 0, this means all observations are censored:
# C1 and C3, give a reasonable CI (-Inf, Inf) or (0, Inf)
# C2,  C4 and C7 give (0, 0)
# C5 and C6 fail

library(ggplot2)

is_valid_ci <- function(ci) {
  !is.null(ci) && length(ci) == 2 && !anyNA(ci)
}

n_seq <- 4:50

coverage_prob_1 <- numeric(50)
coverage_prob_2 <- numeric(50)
coverage_prob_3 <- numeric(50)
coverage_prob_4 <- numeric(50)
coverage_prob_5 <- numeric(50)
coverage_prob_6 <- numeric(50)
coverage_prob_7 <- numeric(50)

for (n in n_seq) {

  n_sim <- 1000
  theta <- 1

  is_in_C1 <- 0
  is_in_C2 <- 0
  is_in_C3 <- 0
  is_in_C4 <- 0
  is_in_C5 <- 0
  is_in_C6 <- 0
  is_in_C7 <- 0

  n_sim_1 <- n_sim
  n_sim_2 <- n_sim
  n_sim_3 <- n_sim
  n_sim_4 <- n_sim
  n_sim_5 <- n_sim
  n_sim_6 <- n_sim
  n_sim_7 <- n_sim
  
  for (i in 1:n_sim) {

    t <- rexp(n = n, rate = 1 / theta)

    C1 <- tryCatch(
      C1_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )
    C2 <- tryCatch(
      C2_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )
    C3 <- tryCatch(
      C3_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )
    C4 <- tryCatch(
      C4_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )
    C5 <- tryCatch(
      C5_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )
    C6 <- tryCatch(
      C6_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )
    C7 <- tryCatch(
      C7_Sundberg(t = t, t_c = 1, alpha = 0.1)$CI,
      error = function(e) NULL
    )

    if (!is_valid_ci(C1)) {
      n_sim_1 <- n_sim_1 - 1
    } else if (theta > C1[1] && theta < C1[2]) {
      is_in_C1 <- is_in_C1 + 1
    }

    if (!is_valid_ci(C2)) {
      n_sim_2 <- n_sim_2 - 1
    } else if (theta > C2[1] && theta < C2[2]) {
      is_in_C2 <- is_in_C2 + 1
    }

    if (!is_valid_ci(C3)) {
      n_sim_3 <- n_sim_3 - 1
    } else if (theta > C3[1] && theta < C3[2]) {
      is_in_C3 <- is_in_C3 + 1
    }

    if (!is_valid_ci(C4)) {
      n_sim_4 <- n_sim_4 - 1
    } else if (theta > C4[1] && theta < C4[2]) {
      is_in_C4 <- is_in_C4 + 1
    }

    if (!is_valid_ci(C5)) {
      n_sim_5 <- n_sim_5 - 1
    } else if (theta > C5[1] && theta < C5[2]) {
      is_in_C5 <- is_in_C5 + 1
    }

    if (!is_valid_ci(C6)) {
      n_sim_6 <- n_sim_6 - 1
    } else if (theta > C6[1] && theta < C6[2]) {
      is_in_C6 <- is_in_C6 + 1
    }

    if (!is_valid_ci(C7)) {
      n_sim_7 <- n_sim_7 - 1
    } else if (theta > C7[1] && theta < C7[2]) {
      is_in_C7 <- is_in_C7 + 1
    }
  }

  coverage_prob_1[n] <- is_in_C1 / n_sim_1
  coverage_prob_2[n] <- is_in_C2 / n_sim_2
  coverage_prob_3[n] <- is_in_C3 / n_sim_3
  coverage_prob_4[n] <- is_in_C4 / n_sim_4
  coverage_prob_5[n] <- is_in_C5 / n_sim_5
  coverage_prob_6[n] <- is_in_C6 / n_sim_6
  coverage_prob_7[n] <- is_in_C7 / n_sim_7
}





results_df <- data.frame(
  n = rep(n_seq, times = 7),
  coverage = c(
    coverage_prob_1[n_seq],
    coverage_prob_2[n_seq],
    coverage_prob_3[n_seq],
    coverage_prob_4[n_seq],
    coverage_prob_5[n_seq],
    coverage_prob_6[n_seq],
    coverage_prob_7[n_seq]
  ),
  method = rep(paste0("C", 1:7), each = length(n_seq))
)


ggplot(results_df, aes(x = n, y = coverage, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40") +
  labs(
    x = "n (sample size)",
    y = "Coverage probability",
    title = "Coverage probability of C1-C7 vs. n",
    color = "Method"
  ) +
  theme_minimal()
