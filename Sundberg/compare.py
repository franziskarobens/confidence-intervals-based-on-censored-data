"""
Compare Ci (Type I specific) and Ci (Sundberg) for type I censored exponential data.

- plot n vs P(non coverage) for both methods (simulation)
- plot n vs mean width of the interval
"""

import numpy as np
import matplotlib.pyplot as plt

import sundberg
from sundberg_type_I import Limits
import sundberg_type_I

method_nr = 5

alpha = 0.1
theta = 1
nu = 0.5
c = nu * theta

n_range = list(range(4, 51))
n_sim = 100_000

rng = np.random.default_rng(2026)

noncov_I, noncov_II = [], []
width_I, width_II = [], []

for n in n_range:
    print(n, end=" ")
    miss_I = miss_II = 0
    w_I, w_II = [], []

    for _ in range(n_sim):
        X = rng.exponential(scale=1 / theta, size=n) # rate = theta
        TTT, R, MLE = sundberg_type_I.get_MLE(X, c, theta_is_rate=True)
        if R == 0:                    # CI undefined
            continue

        C_func = getattr(sundberg_type_I, f"C{method_nr}")
        C_I = C_func(n, c, X, alpha)
        
        c_lower, c_upper = sundberg.CI(method_nr = method_nr, r=R, alpha=alpha, theta_is_rate=True)
        C_II = Limits(c_lower * MLE, c_upper * MLE)
        
        if not (C_I.lower <= theta <= C_I.upper):
            miss_I += 1
        if not (C_II.lower <= theta <= C_II.upper):
            miss_II += 1

        w_I.append(C_I.upper - C_I.lower)
        w_II.append(C_II.upper - C_II.lower)

    noncov_I.append(miss_I / n_sim)
    noncov_II.append(miss_II / n_sim)
    width_I.append(np.mean(w_I))
    width_II.append(np.mean(w_II))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))

ax[0].plot(n_range, noncov_I, color="red", label=f"C{method_nr} Type I")
ax[0].plot(n_range, noncov_II, color="blue", label=f"C{method_nr} Sundberg")
ax[0].axhline(alpha, color="gray", ls="--", label=r"nominal $\alpha$")
ax[0].set_xlabel("n")
ax[0].set_ylabel("P(non coverage)")
ax[0].legend()

ax[1].plot(n_range, width_I, color="red", label=f"C{method_nr} Type I")
ax[1].plot(n_range, width_II, color="blue", label=f"C{method_nr} Sundberg")
ax[1].set_xlabel("n")
ax[1].set_ylabel("mean width")
ax[1].legend()

plt.tight_layout()
plt.show()
