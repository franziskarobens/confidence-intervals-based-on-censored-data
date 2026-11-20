"""
C6 is exact for type II censored data and only asymptotic for type I.

Compare, for the rate theta of an exponential distribution:
  (I)  type I censoring at time c, C6 with the random number of failures R
  (II) type II censoring after r = floor(n/2) failures, C6 with fixed r

- plot n vs P(non coverage) for both settings
- plot n vs mean width of the interval
"""

from functools import lru_cache

import numpy as np
import matplotlib.pyplot as plt

import sundberg

# ---------------------------------------------------------------- parameters
alpha = 0.1
theta = 1.0                      # rate
nu = 0.5
c = nu * theta                   # censoring time (type I)

n_range = list(range(4, 51))
n_sim = 100_000
rng = np.random.default_rng(2026)


# ------------------------------------------------------------------- helpers
@lru_cache(maxsize=None)
def c_constants(r: int, method_nr=6):
    """Quotient limits (c_lower, c_upper) of C6; they depend only on r and alpha."""
    C_func = getattr(sundberg, f"C{method_nr}")
    lim = C_func(r=r, alpha=alpha)
    return float(lim.c_lower), float(lim.c_upper)


def simulate_type_I(n):
    """Return (miss, width) arrays for type I data. Runs with R = 0 are dropped."""
    X = rng.exponential(scale=1 / theta, size=(n_sim, n))
    R = (X < c).sum(axis=1)
    TTT = np.minimum(X, c).sum(axis=1)

    valid = R > 0
    R, TTT = R[valid], TTT[valid]
    mle = R / TTT                                    # rate MLE

    lo_tab = np.array([np.nan] + [c_constants(k)[0] for k in range(1, n + 1)])
    hi_tab = np.array([np.nan] + [c_constants(k)[1] for k in range(1, n + 1)])
    lower, upper = lo_tab[R] * mle, hi_tab[R] * mle
    return ~((lower <= theta) & (theta <= upper)), upper - lower


def simulate_type_II(n, r):
    """Return (miss, width) arrays for type II data (stop at the r-th failure)."""
    X = np.sort(rng.exponential(scale=1 / theta, size=(n_sim, n)), axis=1)
    TTT = X[:, :r].sum(axis=1) + (n - r) * X[:, r - 1]
    mle = r / TTT

    c_lo, c_hi = c_constants(r)
    lower, upper = c_lo * mle, c_hi * mle
    return ~((lower <= theta) & (theta <= upper)), upper - lower


# ---------------------------------------------------------------- simulation
noncov_I, noncov_II, width_I, width_II = [], [], [], []

for n in n_range:
    print(n, end=" ", flush=True)

    miss, w = simulate_type_I(n)
    noncov_I.append(miss.mean())
    width_I.append(w.mean())

    miss, w = simulate_type_II(n, r=n // 2)
    noncov_II.append(miss.mean())
    width_II.append(w.mean())
print()

# -------------------------------------------------------------------- plot

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))

ax[0].plot(n_range, noncov_I, color="red", label="type I (R random)")
ax[0].plot(n_range, noncov_II, color="blue", label=r"type II ($r=\lfloor n/2\rfloor$)")
ax[0].axhline(alpha, color="gray", ls="--", label=r"nominal $\alpha$")
ax[0].set_xlabel("n")
ax[0].set_ylabel("P(non coverage)")
ax[0].legend()

ax[1].plot(n_range, width_I, color="red", label="type I (R random)")
ax[1].plot(n_range, width_II, color="blue", label=r"type II ($r=\lfloor n/2\rfloor$)")
ax[1].set_xlabel("n")
ax[1].set_ylabel("mean width")
ax[1].legend()

plt.tight_layout()
plt.show()