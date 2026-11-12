

theta <- 2
c <- 0.4
r <- 5
n <- 100

X <- rexp(n = n, rate = theta)
X_sorted <- sort(X)

# Type I and Type II censoring
X_I_censored  <- pmin(X, c)
delta_I <- X <= c

X_II_censored <- c(X_sorted[1:r], rep(X_sorted[r], n-r))
delta_II <- X <= X_sorted[r]

# Type-I hybrid: stop at min(X_(r), c)
tau_hybrid_I   <- min(X_sorted[r], c)     # min(X_(r), c)
X_hybrid_I     <- pmin(X, tau_hybrid_I)   # censored units are recorded at the stopping time
delta_hybrid_I <- X <= tau_hybrid_I       # failure indicators

# Type-II hybrid: stop at max(X_(r), c)
tau_hybrid_II   <- max(X_sorted[r], c)
X_hybrid_II     <- pmin(X, tau_hybrid_II)
delta_hybrid_II <- X <= tau_hybrid_II

# MLE = number of failures / total time on test
mle <- function(x_obs, delta) sum(delta) / sum(x_obs)

n / sum(X)
mle(X_I_censored, delta_I)
mle(X_II_censored, delta_II)
mle(X_hybrid_I,  delta_hybrid_I)
mle(X_hybrid_II, delta_hybrid_II)




## PIVOTING THE CDF

alpha <- 0.1
theta_true <- 1
n_sim <- 100000
X <- rexp(n = n_sim, rate = theta_true)

CI <- matrix(nrow = n_sim, ncol = 2)
in_CI <- 0

for (i in 1:n_sim) {
  lower_limit <- - log(1 - alpha / 2) / X[i]
  upper_limit <- - log(alpha / 2) / X[i]
  CI[i, ] <- c(lower_limit, upper_limit)
  if (theta_true >= lower_limit & theta_true <= upper_limit) {
    in_CI <- in_CI + 1
  }
}

in_CI / n_sim


# for type II censored 

n <- 10
n_sim <- 10000
r <- 4 # r < n

in_CI_type_II <- 0

for (i in 1:n_sim) {
  X <- rexp(n = n, rate = theta_true)
  sorted_X <- sort(X)
  T <- sum(sorted_X[1:r]) + (n-r) * sorted_X[r]

  g <- qgamma(p = c(alpha/2, 1 - alpha/2), shape = r, rate = 1)
  CI <- c(g[1] / T, g[2] / T)

  if (theta_true >= CI[1] & theta_true <= CI[2]) {
    in_CI_type_II <- in_CI_type_II + 1
  }
}

in_CI_type_II / n_sim


# for type I censored 

c <- 0.5
theta_true <- 1
n <- 10
n_sim <- 1000
alpha <- 0.1


dens <- function(x, k, n, cc, theta) {
  sapply(x, function(xx) {
    j <- 0:min(k, floor(xx / cc))
    s <- sum((-1)^j * choose(k, j) * (xx - j * cc)^(k - 1))
    choose(n, k) * exp(-theta * cc * (n - k)) * theta^k * exp(-theta * xx) /
      factorial(k - 1) * s
  })
}

# P(T <= t | R >= 1), T = R/S
F_cond <- function(t, theta, n, cc) {
  tot <- 0
  for (k in 1:n) {
    L <- max(0, k / t - (n - k) * cc)
    U <- k * cc
    if (L >= U) next
    kinks <- cc * (1:k)
    br <- sort(unique(c(L, kinks[kinks > L & kinks < U], U)))  # split at kinks of the density
    for (i in seq_len(length(br) - 1))
      tot <- tot + integrate(dens, br[i], br[i + 1],
                             k = k, n = n, cc = cc, theta = theta)$value
  }
  tot / (1 - exp(-n * theta * cc))
}

ci_type1 <- function(R, S, n, cc, alpha = 0.05) {
  t <- R / S
  f <- function(theta, p) F_cond(t, theta, n, cc) - p
  lo <- tryCatch(uniroot(f, c(1e-4, 100), p = 1 - alpha / 2)$root,
                 error = function(e) 0)      # no solution: lower limit is 0
  hi <- tryCatch(uniroot(f, c(1e-4, 100), p = alpha / 2)$root,
                 error = function(e) Inf)    # no solution: upper limit is Inf
  c(lower = lo, upper = hi)
}

in_CI_type_I <- 0
n_valid <- 0            # number of simulations with R >= 1

for (i in 1:n_sim) {
  X <- rexp(n = n, rate = theta_true)
  sorted_X <- sort(X)
  R <- sum(X <= c)
  T <- sum(sorted_X[1:R]) + (n-R) * c

  if (R == 0) next                      # the interval is only defined given R >= 1

  S <- sum(pmin(X, c))                 # total time on test
  n_valid <- n_valid + 1

  CI <- ci_type1(R = R, S = S, n = n, cc = c, alpha = alpha)

  if (theta_true >= CI[1] & theta_true <= CI[2]) {
    in_CI_type_I <- in_CI_type_I + 1
  }
}

in_CI_type_I / n_sim



# ------------- Paper von Hannig (Figure 1) ----------------- 

theta <- 100
y <- runif(n = 25, theta, theta^2)                      # your sample
n <- length(y)
theta <- seq(98, 110, length.out = 5000)
lo <- sqrt(max(y)); hi <- min(y)      # support of the GFD

# unnormalized density (here with J = sum of |dy_i/dtheta|)
dens <- function(t) {
  J <- sapply(t, function(th)
    sum(abs(1 + (2*th - 1) * (y - th) / (th^2 - th))))
  (t^2 - t)^(-n) * J
}

dens_vals <- ifelse(theta > lo & theta < hi, dens(theta), 0)
R  <- cumsum(dens_vals) / sum(dens_vals)   # fiducial cdf on the grid
cc <- abs(1 - 2 * R)

plot(theta, cc, type = "l", ylim = c(0, 1),
     xlab = expression(theta), ylab = "confidence level")
abline(h = 0.95, lty = 3)  





# Woher kommt die Oszillation?

# non cov probs for nu = 0.5, alpha = 0.1, C1 Methode und theta = 1, Digits = 1
x <- c(0.0, 0.17313580830980557, 0.20394121684215483, 0.14139763316243237, 0.10295028505391324, 0.20704519536340515, 0.1724742945570752, 0.12452503998148201, 0.14351749155497906, 0.17800069513781624, 0.13401045462897432, 0.12016853867097921, 0.1627667031679717, 0.13749991045725818, 0.1143988634608436, 0.1433438900768551, 0.1374790634212555, 0.11295987019802334, 0.12922659638054346, 0.1341689580180192, 0.11229845150110085, 0.12021918891413463, 0.12913498857878902, 0.11149625710576228, 0.11445758195992395, 0.12384894049288386, 0.11030503904224731, 0.1105708566925031, 0.11902168545873149, 0.10875578918629769, 0.10775586880201304, 0.11484781599303107, 0.10698570467396129, 0.1055643414330576, 0.11129582735066212, 0.10513368454036731, 0.10374806515131707, 0.1082649643603879, 0.10330174136581038, 0.1021693839995933, 0.10565013764512776, 0.10155253655999175, 0.10075091996463015, 0.10336298452846215, 0.0999188019892977, 0.09944759506568746, 0.10133547464113568)
6 9 13 
n <- seq(4, 50, length.out = 1000)
exp_R <- n * (1 - exp(-0.5))


library(reticulate)

sundberg_dir <- "C:/Users/Franzi/Documents/mathe_studium/masterarbeit/confidence-intervals-based-on-censored-data/Sundberg"

sundberg <- import_from_path("sundberg", path = sundberg_dir)

c_lower <- numeric(length(n))
c_upper <- numeric(length(n))


for (i in seq_along(n)) {
  res <- sundberg$CI(method_nr = 1L, r = n[i], alpha = 0.1, theta_is_rate = FALSE)
  c_lower[i] <- res[[1]]
  c_upper[i] <- res[[2]]
}

plot(n, x, type = "l")
plot(n, c_lower)
plot(n, c_upper)

