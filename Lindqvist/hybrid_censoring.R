# ---- Plot eq. (8) ----------------------------
# DGF of MLE for hybrid type II censoring

n <- 8
true_theta <- 1
c <- 1
r <- 3

X <- rexp(n = n, rate = true_theta)
sort_X <- sort(X)
R_obs <- sum(X <= c)
if (R_obs < r) T_obs <- r / (sum(sort_X[1:r]) + (n-r) * sort_X[r])
if (R_obs >= r) T_obs <- R_obs / (sum(sort_X[1:R_obs]) + (n-R_obs) * c)


U <- rexp(n = n, rate = 1)

tau_8 <- function(U, r, c) {
  theta_grid <- seq(0, max(U / c), length.out = 1000)
  tau <- numeric(1000)
  U_sort <- sort(U)

  for (i in 1:1000) {
    theta <- theta_grid[i]
    R <- sum(U / theta <= c)

    if (R < r){
      tau[i] <- r * theta / (sum(U_sort[1:r]) + (n-r) * U_sort[r])
    }
    else {
      tau[i] <- R / (sum(U_sort[1:r]) / theta + (n-R) * c)
    }
  }

  return (list(theta = theta_grid, 
  tau = tau))
}

result <- tau_8(U, r, c)

plot(x = result$theta, y = result$tau)
abline(v = U)
abline(h = T_obs, col = "red")



# DGF of MLE for hybrid type I censoring


n <- 10
true_theta <- 1
c <- 1
r <- 5

X <- rexp(n = n, rate = true_theta)
sort_X <- sort(X)
R_obs <- sum(X <= c)
if (R_obs >= r) T_obs <- r / (sum(sort_X[1:r]) + (n-r) * sort_X[r])
if (R_obs < r) T_obs <- R_obs / (sum(sort_X[1:R_obs]) + (n-R_obs) * c)



U <- rexp(n = n, rate = 1)

tau_8 <- function(U, r, c) {
  theta_grid <- seq(0, max(U / c), length.out = 1000)
  tau <- numeric(1000)
  U_sort <- sort(U)

  for (i in 1:1000) {
    theta <- theta_grid[i]
    R <- sum(U / theta <= c)

    if (R >= r){
      tau[i] <- r * theta / (sum(U_sort[1:r]) + (n-r) * U_sort[r])
    }
    else {
      tau[i] <- R / (sum(U_sort[1:r]) / theta + (n-R) * c)
    }
  }

  return (list(theta = theta_grid, 
  tau = tau))
}

result <- tau_8(U, r, c)

plot(x = result$theta, y = result$tau)
abline(v = U)
abline(h = T_obs, col = "red")

