# Description
# TBD


# returns MLE given by Sundberg in equation (2) in chapter 2.2
# @ params
# t = c(t_1, ..., t_n) with t_i ~ Exp
# t_c censoring time
mle_Sundberg <- function(t, t_c) {
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  t_obs <- pmin(t, t_c)

  return(sum(t_obs) / N)
}

# returns CI regarding to Sundberg (see C1)
# @params
# t is a n dimensional vector of observations whereas we assume 
# they are exponentially distributed
# t_c censoring time 
# alpha confidence niveau for the CI

C1_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  # maximum likelihood estimations
  theta_hat <- mle_Sundberg(t, t_c)
  
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  z <- qnorm(p = 1 - alpha / 2)
  
  ci <- c(theta_hat * (1 - z / sqrt(N)), theta_hat * (1 + z / sqrt(N)))
  
  return(unname(ci))
}


C2_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  # maximum likelihood estimations
  theta_hat <- mle_Sundberg(t, t_c)
  
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  z <- qnorm(p = 1 - alpha / 2)
  
  ci <- c(theta_hat / (1 + z / sqrt(N)), theta_hat / (1 - z / sqrt(N)))
  
  return(unname(ci))
}

C3_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  # maximum likelihood estimations
  theta_hat <- mle_Sundberg(t, t_c)
  
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  z <- qnorm(p = 1 - alpha / 2)
  
  ci <- c(theta_hat * exp(- z / sqrt(N)), theta_hat * exp(z / sqrt(N)))
  
  return(unname(ci))
}

C4_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  # maximum likelihood estimations
  theta_hat <- mle_Sundberg(t, t_c)
  
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  z <- qnorm(p = 1 - alpha / 2)
  
  ci <- c(theta_hat / (1 + z / (3*sqrt(N)))^3, theta_hat / (1 - z / (3*sqrt(N)))^3)
  
  return(unname(ci))
}


C5_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  theta_hat <- mle_Sundberg(t, t_c)
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
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
  
  ci <- c(b1, b2)
  
  return(unname(ci))
}



C6_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  theta_hat <- mle_Sundberg(t, t_c)
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  q1 <- qchisq(p = 1 - alpha / 2, df = 2 * N)
  q2 <- qchisq(p = alpha / 2, df = 2 * N)
  
  ci <- c(2 * N * theta_hat / q1, 2 * N * theta_hat / q2)
  
  return(unname(ci))
}

C7_Sundberg <- function(t, t_c, alpha = 0.1) {
  
  theta_hat <- mle_Sundberg(t, t_c)
  N <- sum(get_delta(data = t, t_c = t_c)) # number of uncensored
  q1 <- qchisq(p = 1 - alpha / 2, df = 2 * N + 1)
  q2 <- qchisq(p = alpha / 2, df = 2 * N + 1)
  
  ci <- c(2 * N * theta_hat / q1, 2 * N * theta_hat / q2)
  
  return(unname(ci))
}


C <- 0
C1_Sundberg(t = ball_bearing_data, t_c = C)
C2_Sundberg(t = ball_bearing_data, t_c = C)
C3_Sundberg(t = ball_bearing_data, t_c = C)
C4_Sundberg(t = ball_bearing_data, t_c = C)
C5_Sundberg(t = ball_bearing_data, t_c = C)
C6_Sundberg(t = ball_bearing_data, t_c = C)
C7_Sundberg(t = ball_bearing_data, t_c = C)

# Observation: confidence intervals C1 and C2 differ a lot from the remaining intervals
# C3 - C7 are relatively similar

# for edge cases:
# C = 0, this means all observations are censored:
# C1 and C3, give a reasonable CI (-Inf, Inf) or (0, Inf)
# C2,  C4 and C7 give (0, 0)
# C5 and C6 fail