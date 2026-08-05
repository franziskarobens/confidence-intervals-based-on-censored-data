
source("plot_ci_limits.R")

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

plot_ci_limits(t = t, alpha = 0.1)
