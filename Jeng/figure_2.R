# Replicate Figure 2 from Jeng

p_f_vals <- c(0.01, 0.1, 0.3, 0.5, 1)

# first NORM procedure

E_r_vals <- 2:100

ci_NORM_Jeng()



coverage <- function(ci_fun, n, p_f, alpha = 0.1, M = 2000) {
  t_c <- log(- log(1 - p_f)) # solution of P(T <= t_c) = p_f, T ~ SEV(0, 1)
  hit <- logical(M)          # FALSE if true parameter not in CI, TRUE otherwise
  ok <- logical(M)           # FALSE if r <= 1, and TRUE if r > 1
  for (m in 1:M) {
    x <- rsev(n)
    r <- sum(x <= t_c)                   # Number of failures (Type I)
    T <- sum(pmin(x, t_c))               # Total time on test
    if (r <= 1) { ok[m] <- FALSE; next } # Skip this case
    #### ci <- ci_fun(T, n, r, t_c, alpha)
    # z. B. NORM_Jeng: blablab

    hit[m] <- ci["lower"] <= 0 && 0 <= ci["upper"]
    ok[m] <- TRUE
  }
  c(coverage = mean(hit[ok]), mc_se = sqrt(mean(hit[ok]) * (1 - mean(hit[ok])) / sum(ok)),
    n_valid = sum(ok))
}

n = E_r / p_f
coverage(ci_fun = NORM_CI, n = n, p_f = 0.1)
NORM_CI(t, t_c, alpha)