# constants c_lower and c_upper form C1, ... , C7 from Sundberg
# assuming parametrization is rate = theta

C1_Sundberg <- function(r, alpha) {
  z <- qnorm(p = 1 - alpha / 2)

  c(c_lower = 1 - z / sqrt(r), 
    c_upper = 1 + z / sqrt(r))
}


C2_Sundberg <- function(r, alpha) {
  if (r < qnorm(p = 1 - alpha / 2)^2) {
    warning("C2 limits aren't valid. need higher r.")
  }

  z <- qnorm(p = 1 - alpha / 2)

  c(c_lower = 1 / (1 + z / sqrt(r)), 
    c_upper = 1 / (1 - z / sqrt(r)))
}


C3_Sundberg <- function(r, alpha) {
  z <- qnorm(p = 1 - alpha / 2)

  c(c_lower = exp(- z / sqrt(r)), 
    c_upper = exp(z / sqrt(r)))
}


C4_Sundberg <- function(r, alpha) {
  z <- qnorm(p = 1 - alpha / 2)
  c(c_lower = 1 / (1 + z / (3 * sqrt(r)))^3, 
    c_upper = 1 / (1 - z / (3 * sqrt(r)))^3)
}


C5_Sundberg <- function(r, alpha) {
  z <- qnorm(1 - alpha / 2)
  t <- z^2 / (2 * r)

  arg <- -exp(-t - 1)          # lays in (-1/e, 0) for t > 0
  c1  <- -lambertW0(arg)       # solution in (0, 1)
  c2  <- -lambertWm1(arg)      # solution in (1, Inf)

  c(c_lower = c1, c_upper = c2)
}


C6_Sundberg <- function(r, alpha) {
  c1 <- qchisq(p = alpha / 2, df = 2 * r)
  c2 <- qchisq(p = 1 - alpha / 2, df = 2 * r)
  
  c(c_lower = c1 / (2 * r), 
    c_upper = c2 / (2 * r))
}


C7_Sundberg <- function(r, alpha) {
  c1 <- qchisq(p = alpha / 2, df = 2 * r + 1)
  c2 <- qchisq(p = 1 - alpha / 2, df = 2 * r + 1)
  
  c(c_lower = c1 / (2 * r), 
    c_upper = c2 / (2 * r))
}


Sundberg_limits <- function(method_nr, r, alpha) {
  method_name <- paste0("C", method_nr, "_Sundberg")
  f <- match.fun(method_name)
  return(f(r, alpha))
}
