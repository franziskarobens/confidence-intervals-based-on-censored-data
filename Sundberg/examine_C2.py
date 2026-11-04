# Discover the problem that for C2 in Sundberg, for small n 
# the interval bounds don't fulfil c_l < c_u
# Plot sample size n vs CI width

import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from sundberg import C2


n = np.arange(0, 11)                    # 0, ...., 10

alpha = 0.01
color = "red" if alpha == 0.1 else "blue"

# Multipliers for C2 from Sundberg (uncensored case: r = n)
with warnings.catch_warnings(), np.errstate(divide="ignore", invalid="ignore"):
    warnings.simplefilter("ignore")    # C2 warns for r < z^2, which is expected here
    c_lower, c_upper = C2(n, alpha)

widths = c_upper - c_lower

# For all n < z^2, it holds c_l > c_u.
z2 = norm.ppf(1 - alpha / 2) ** 2

fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(n, widths, "o", color=color, markersize=5)

# add vertical line where width becomes positive
ax.axvline(z2, color=color)

ax.set_xticks(n)
ax.set_xlabel(r"$n$")
ax.set_ylabel(r"$c_u - c_l$")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(alpha=0.3)

fig.tight_layout()
plt.show()