"""Use Theorem 1 and Corollary 1, 2 to create a CI for 
exponentially distributed data."""

import numpy as np
from scipy import stats


# CI for uncensored exponential Exp(1/theta) data
def CI_corollary_2(theta, T, n, alpha, is_theta_rate=True):
    b = T / stats.gamma.ppf(alpha / 2, a=n, scale=1)
    a = T / stats.gamma.ppf(1 - alpha / 2, a=n, scale=1)
    
    if is_theta_rate:
        return a, b
    return 1 / b, 1 / a


alpha = 0.1
n = 100
n_sim = 10000
theta = 1

covered = 0

for i in range(n_sim):
    rng = np.random.default_rng(i)
    X = rng.exponential(scale=1/theta, size=n)     # Exp(rate = theta)
    T = np.sum(X)

    lower, upper = CI_corollary_2(theta, T = T, n = n, alpha = alpha, is_theta_rate = True)
    print(f"n_sim = {i}, lower = {lower}, higher = {upper}, T = {T}")

    covered += lower <= theta <= upper


print(covered / n_sim)   # should be close to 1 - alpha



"""One cannot find a and b for censored data, because to get a and b one need the 
distribution of bar theta and undeline theta."""
# TODO: find a possibility to get bar theta and underline theta
# use the DGF and find inverse (should be possible bc we have DGF given)