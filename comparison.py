"""
Vergleich von Konfidenzintervall-Methoden für Type-I-zensierte
exponentialverteilte Daten: Coverage Probability, mittlere Breite, Laufzeit.

Modell: X_i = U_i / theta, U_i ~ Exp(1), Zensierung bei c.
Statistik T = (#Ausfälle) / (Gesamttestzeit)  (= ML-Schätzer von theta).
"""

import time
from dataclasses import dataclass
from typing import Callable
from Sundberg import sundberg
from Lindqvist import lindqvist_monte_carlo

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# create data
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Sample:
    n: int          # sample size
    c: float        # censoring time
    U: float        # random vector on which the observation is based
    X: float        # observation
    R: int          # number of observed values
    TTT: float      # total time on test
    MLE: float      # MLE = R / TTT (for Exp(rate = 1/theta))


def draw_sample(theta, n, nu, rng, theta_is_rate=True) -> Sample:
    """
    Draw a Type I censored exponential sample.
    theta_is_rate=False: X ~ Exp(rate = 1/theta), i.e. X = theta * U
    theta_is_rate=True:  X ~ Exp(rate = theta), i.e. X = U / theta
    """
    scale = 1 / theta if theta_is_rate else theta   # mean of X
    c = nu * scale                                  # censoring time
    U = rng.exponential(scale=1.0, size=n)          # U ~ Exp(1)
    X = scale * U
    R = int(np.sum(X < c))
    TTT = float(np.sum(np.minimum(X, c)))

    # MLE of theta; undefined if no failure is observed
    if R == 0:
        MLE = np.nan
    else:
        MLE = R / TTT if theta_is_rate else TTT / R

    return Sample(n, c, U, X, R, TTT, MLE)

# --- Define all methods observed -----------------------------------------

# every method is (int: r, float: alpha) -> (float: c_l)
Method = Callable[[int, float], tuple[float, float]]

def make_sundberg(method_nr: int) -> Method:
    def method(sample, alpha):
        c_lower, c_upper = sundberg.CI(method_nr=method_nr, r=sample.R, alpha=alpha, theta_is_rate=True)
        return c_lower * sample.MLE, c_upper * sample.MLE 
    return method

def make_lindqvist_monte_carlo(m: int, k: int) -> Method:
    def method(sample, alpha):
        level = 2 * k / (m + 1)
        if alpha != level:
            raise ValueError(f"This method has level {level}, got alpha={alpha}.")
        return lindqvist_monte_carlo.CI_corollary_4(T=sample.MLE, n=sample.n, c=sample.c, k=k, m=m)
    return method

# (m, k) pairs, alpha = 2k/(m+1)
DGF_CONFIGS = [(199, 1), (999, 5), (99, 5), (999, 50)]   # alpha = 0.01, 0.01, 0.1, 0.1

# define all methods
METHODS: dict[str, Method] = {
    **{f"C{i}_Sundberg": make_sundberg(i) for i in range(1, 8)},
    **{f"Lindqvist_MC_k{k}_m{m}": make_lindqvist_monte_carlo(m, k) for m, k in DGF_CONFIGS},
}

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def run_scenario(methods, theta, n, nu, alpha, n_sim, seed):
    rng = np.random.default_rng(seed)
    stats_ = {name: dict(n_ok=0, n_failed=0, cover=0, miss_low=0, miss_up=0,
                         width=0.0, time=0.0)
              for name in methods}
    valid = 0
    while valid < n_sim:
        sample = draw_sample(theta, n, nu, rng)
        if sample.R == 0:          # T = 0: CI not defined -> condition on r > 0
            continue
        valid += 1
        for name, method in methods.items():
            d = stats_[name]
            t0 = time.perf_counter()
            try:
                lower, upper = method(sample, alpha)
            except ValueError:     # e.g. alpha not supported by this method
                continue
            d["time"] += time.perf_counter() - t0
            d["n_ok"] += 1
            d["cover"] += lower <= theta <= upper
            d["miss_low"] += theta < lower      # theta under CI
            d["miss_up"] += theta > upper       # theta over CI
            d["width"] += upper - lower

    rows = []
    for name, d in stats_.items():
        if d["n_ok"] == 0:         # method never produced a CI: no row
            continue
        p = d["cover"] / n_sim
        rows.append({
            "method": name,
            "theta": theta,
            "n": n,
            "nu": nu,
            "c": nu * theta,
            "P(censored)": np.exp(-nu),
            "coverage": p,
            "mc_se": np.sqrt(p * (1 - p) / n_sim),
            "miss_below": d["miss_low"] / n_sim,
            "miss_above": d["miss_up"] / n_sim,
            "mean_width": d["width"] / n_sim,
            "ms_per_ci": 1000 * d["time"] / n_sim,
        })
    return rows


def main():
    alpha = 0.1
    theta = 1.0
    n_sim = 1000

    rows = []
    for n in (10, 50, 100, 1000): # 10, 50, 100, 1000
        for nu in (0.5, 1, 2): # 0.5, 1, 2
            print(f"Scenario n={n}, nu={nu} ...", flush=True)
            rows += run_scenario(METHODS, theta, n, nu, alpha, n_sim,
                                 seed=1000 * n + int(10 * nu))

    df = pd.DataFrame(rows)
    df.to_csv(f"data_comparison/ci_comparison_alpha{alpha}.csv", index=False)


if __name__ == "__main__":
    main()
