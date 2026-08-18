
g_log <- function(x) {
  g(x) <- log(x)
  g_inv(x) <- exp(x)
  g_der(x) <- 1/x
  
  # TODO return all functions
}

ci_tnorm <- function(t, t_c, g, alpha) {
  z <- qnorm(q = 1 - alpha/2)
  
}