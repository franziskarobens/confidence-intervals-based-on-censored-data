# Replicate the Figures from Sundberg


source("Sundberg/plot_CI_limits.R")

theta <- 100
lambda <- 1 / theta
t <- rexp(n = 10000, lambda)

C <- 0.5 * theta
C <- theta

C1_Sundberg(t, t_c = C)
C2_Sundberg(t, t_c = C)
C3_Sundberg(t, t_c = C)
C4_Sundberg(t, t_c = C)
C5_Sundberg(t, t_c = C)
C6_Sundberg(t, t_c = C)
C7_Sundberg(t, t_c = C)


t <- rsev(n = 100, mu = 4, sigma = 3)
C <- 4

ci_NORM_Jeng(t, t_c = C) # hier wird die fisher nicht berechnet weil die hesse Nullmatrix ist 
# TBD das fixen

plot_CI_limits(t = t, alpha = 0.1)
