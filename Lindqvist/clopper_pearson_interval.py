import numpy as np
from scipy.stats import beta

# --- Example: Clopper-Pearson Interval (Section 2.1 in Lindqvist) --------------

n = 100
alpha = 0.05
theta = 0.5
n_simulations = 10**7

rng = np.random.default_rng()

# a[t] = upper limit for T = t, b[t] = lower limit for T = t  (t = 0, ..., n)
a = np.empty(n + 1)
b = np.empty(n + 1)

a[n] = 1.0   # a(n) = 1
b[0] = 0.0   # b(0) = 0

for t in range(0, n):          # t = 0, ..., n-1
    a[t] = beta.ppf(1 - alpha, t + 1, n - t)

for t in range(1, n + 1):      # t = 1, ..., n
    b[t] = beta.ppf(alpha, t, n + 1 - t)

T = rng.binomial(n=n, p=theta, size=n_simulations)

# upper CI: [0, a[T]]
# lower CI: [b[T], 1]
theta_in_upper_CI = theta <= a[T]
theta_in_lower_CI = theta >= b[T]

# should be >= 1 - alpha (Clopper-Pearson is conservative)
print(theta_in_upper_CI.mean())
print(theta_in_lower_CI.mean())