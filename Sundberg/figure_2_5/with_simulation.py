"""Replicate Figures 2 - 5 from Sundberg using Monte Carlo simulation
and compare with the theoretical values from the CSV files."""
# Observation: the oszillation is also part of the simulation. So it seems to be correct
# that the non-coverage prob is not monotone decreasing

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sundberg

n_monte_carlo = 100_000
nu = 0.5
theta = 1
alpha = 0.1            # level of the two-sided CI passed to sundberg.CI
method_nr = 7
Digits = 100           # as in the file name of the theoretical results

rng = np.random.default_rng(1)

n_range = range(4, 51)
non_cov_probs = []
frac_r0 = []

c = nu * theta

for n in n_range:
    lo = np.full(n + 1, np.nan)
    hi = np.full(n + 1, np.nan)
    for r_val in range(1, n + 1):
        lo[r_val], hi[r_val] = sundberg.CI(method_nr=method_nr, r=r_val,
                                           alpha=alpha, theta_is_rate=False)

    x = rng.exponential(scale=theta, size=(n_monte_carlo, n))
    censored = np.minimum(x, c)
    r = (x <= c).sum(axis=1)

    valid = r > 0
    mle = censored[valid].sum(axis=1) / r[valid]

    covered = (lo[r[valid]] * mle <= theta) & (theta <= hi[r[valid]] * mle)

    non_cov_probs.append(1 - covered.mean())
    frac_r0.append(1 - valid.mean())

# ---------------------------------------------------------------------------
# Theoretical results from the CSV (two-sided CI, same parameters)
# ---------------------------------------------------------------------------
csv_dir = Path(__file__).resolve().parent / "non_coverage_probs"
csv_file = csv_dir / f"non_coverage_two_sided_nu{nu}_alpha{alpha}_theta{theta}_Digits{Digits}.csv"

if not csv_file.exists():
    available = "\n  ".join(p.name for p in sorted(csv_dir.glob("*.csv")))
    raise FileNotFoundError(f"{csv_file} not found.\nAvailable files:\n  {available}")

df_theory = pd.read_csv(csv_file)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(n_range, non_cov_probs, marker="o", markersize=3,
        label=f"C{method_nr} (Monte Carlo)")
ax.plot(df_theory["n"], df_theory[f"C{method_nr}"], linewidth=1.5,
        label=f"C{method_nr} (theoretical)")
ax.set_ylim(0.07, 0.2)
ax.axhline(alpha, linestyle="--", color="grey", label=r"nominal $\alpha$")
ax.set_xlabel("n")
ax.set_ylabel("non-coverage probability")
ax.set_title(rf"two-sided CI: $\nu$ = {nu:g}, $\alpha$ = {alpha:g}, $\theta$ = {theta:g}")
ax.legend()
plt.show()