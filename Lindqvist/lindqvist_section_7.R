

# get MLE in hybrid censoring case (see equation (8) in Lindqvist)
get_tau_eq_8 <- function(u, theta, c, r) {

  R <- sum(u < c)
  u_sorted <- sort(u)

  if (R < r) {
    denom <- sum(u_sorted[1:r] / theta) + (n-r) * (u_sorted[r] / theta)
    return (r / denom)
  } else {
    denom <- sum(u_sorted[1:R] / theta) + (n-R) * c
    return(R / denom)
  }
}


n <- 2
r <- 1
u <- rexp(n)
theta <- 2
c <- 1

x_grid <- seq(0, 1, length.out = 1000)
y_grid <- sapply(x_grid, 
  FUN = get_tau_eq_8, 
  u = u, 
  c = c, 
  r = r)


plot(x_grid, y_grid)
abline(v = min(u) / c)
abline(h = 1 / (2 * c)) # paper says that (U(1)/c, 1/2c) is a point on the graph

# Question: why is this plot an exp curve? 
# I was expecting a piecewise constant strictly monotonous graph