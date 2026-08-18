# Density
dtexp <- function(x, rate = 1, a = 0, b = Inf) {
  dexp(x, rate) / (pexp(b, rate) - pexp(a, rate)) * (x >= a & x <= b)
}

# CDF
ptexp <- function(x, rate = 1, a = 0, b = Inf) {
  (pexp(x, rate) - pexp(a, rate)) / (pexp(b, rate) - pexp(a, rate))
}

# Quantile function (inverse CDF)
qtexp <- function(p, rate = 1, a = 0, b = Inf) {
  Fa <- pexp(a, rate)
  Fb <- pexp(b, rate)
  qexp(Fa + p * (Fb - Fa), rate)
}

# Random generation via inverse transform sampling
rtexp <- function(n, rate = 1, a = 0, b = Inf) {
  qtexp(runif(n), rate, a, b)
}
