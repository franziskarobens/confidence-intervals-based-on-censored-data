"""
Monte Carlo study of the asymptotic distribution of the standardized MLE of the
exponential rate parameter lambda under different censoring schemes:

    uncensored, Type I, Type II, Hybrid Type I, Hybrid Type II

Idea of the refactoring
-----------------------
All schemes share the same structure: the experiment is stopped at some
(possibly random) time `stop`. Then

    N = number of failures observed up to `stop`   (= D in the notes below)
    S = total time on test = sum_i min(t_i, stop)
    MLE = N / S

The schemes differ only in how `stop` is chosen (t_(m) = m-th order statistic):

    Type I       : stop = C
    Type II      : stop = t_(m)
    Hybrid Type I : stop = min(t_(m), C)
    Hybrid Type II: stop = max(t_(m), C)

So one vectorised simulation function replaces the four loops of the R script.

Observation
-----------
The Fisher information is I(lambda) = D / lambda^2, where D depends on the
censoring scheme. With the correct D the standardized MLE
sqrt(D) * (lambda_hat - lambda) / lambda_hat is approximately N(0, 1)
(lambda_hat ~ N(lambda, lambda^2 / D)).

    Uncensored    : D = n
    Type I        : D = N (number of observations before C)
    Type II       : D = m (fixed, no randomness)
    Hybrid Type I : D = min(m, N), N = number of observations before C
    Hybrid Type II: D = max(m, N)

So D is always the number of actually observed failures. If we divide by n
instead of D, we still get a normal distribution, but

    sqrt(n) * (lambda_hat - lambda) / lambda ~ N(0, n / D)  with n / D > 1,

i.e. the variance of the standardized MLEs is too large.

Problem: D has to be tracked. Easy for Type II (D = m), but D is random for
all other schemes.

Where is the histogram of N centred?
------------------------------------
Number of failures before C: N_C ~ Binomial(n, 1 - exp(-lambda * C)), i.e. mean
n * (1 - exp(-lambda * C)) = 1000 * (1 - e^-1) ~ 632.1 and sd ~ 15.3 for the
settings below. Hybrid I: N = min(m, N_C) (cut off from above at m).
Hybrid II: N = max(m, N_C) (cut off from below at m).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

try:
    from tqdm import tqdm
except ImportError:  # progress bar is optional
    def tqdm(iterable, **kwargs):
        return iterable


# --------------------------------------------------------------------------
# Scheme definitions
# --------------------------------------------------------------------------

SCHEMES = {
    #  name       label                  needs C  needs m
    "type1":   ("Type I censoring",        True,   False),
    "type2":   ("Type II censoring",       False,  True),
    "hybrid1": ("Hybrid Type I censoring", True,   True),
    "hybrid2": ("Hybrid Type II censoring", True,  True),
}


def stopping_times(t, scheme, C=None, m=None):
    """Stopping time of the experiment for each row of `t` (shape: rows x n)."""
    rows = t.shape[0]
    if scheme == "type1":
        return np.full(rows, C, dtype=float)

    # m-th smallest value per row (np.partition is faster than a full sort)
    t_m = np.partition(t, m - 1, axis=1)[:, m - 1]

    if scheme == "type2":
        return t_m
    if scheme == "hybrid1":
        return np.minimum(t_m, C)
    if scheme == "hybrid2":
        return np.maximum(t_m, C)
    raise ValueError(f"Unknown scheme: {scheme!r}")


def censored_mle(t, stop):
    """Number of observed failures N and MLE N / S for each row of `t`."""
    stop = stop[:, None]
    N = (t <= stop).sum(axis=1)             # observed (uncensored) failures
    S = np.minimum(t, stop).sum(axis=1)     # total time on test
    return N, N / S


# --------------------------------------------------------------------------
# Simulation
# --------------------------------------------------------------------------

@dataclass
class SimulationResult:
    scheme: str
    lam: float
    n: int
    C: float | None
    m: int | None
    uncensored_mles: np.ndarray
    censored_mles: np.ndarray
    censored_N: np.ndarray                      # D for the censored MLE
    stopped_by_time: np.ndarray | None = None   # True if stopped at C (not at t_(m))

    @property
    def label(self) -> str:
        return SCHEMES[self.scheme][0]

    @property
    def uncensored_N(self) -> np.ndarray:
        return np.full(len(self.uncensored_mles), self.n)

    @property
    def parameter_text(self) -> str:
        parts = []
        if self.C is not None:
            parts.append(f"C = {self.C}")
        if self.m is not None:
            parts.append(f"m = {self.m}")
        return ", ".join(parts)


def simulate(scheme, lam, n, n_simulations, C=None, m=None,
             chunk_size=500, rng=None):
    """
    Simulate `n_simulations` samples of n iid Exp(lam) lifetimes and compute the
    uncensored and the censored MLE for each of them.

    Samples are generated in chunks to keep memory small.
    """
    if scheme not in SCHEMES:
        raise ValueError(f"scheme must be one of {list(SCHEMES)}, got {scheme!r}")
    _, needs_C, needs_m = SCHEMES[scheme]
    if needs_C and C is None:
        raise ValueError(f"scheme {scheme!r} needs the parameter C")
    if needs_m:
        if m is None:
            raise ValueError(f"scheme {scheme!r} needs the parameter m")
        if not 1 <= m < n:
            raise ValueError(f"m must satisfy 1 <= m < n. m = {m} and n = {n}")
    C = C if needs_C else None
    m = m if needs_m else None

    rng = np.random.default_rng() if rng is None else rng

    uncensored = np.empty(n_simulations)
    censored = np.empty(n_simulations)
    N_c = np.empty(n_simulations, dtype=int)
    by_time = np.empty(n_simulations, dtype=bool) if C is not None else None
    n_without_time_censoring = 0

    for start in tqdm(range(0, n_simulations, chunk_size), desc=scheme):
        rows = min(chunk_size, n_simulations - start)
        block = slice(start, start + rows)

        t = rng.exponential(scale=1 / lam, size=(rows, n))   # rexp(rate = lam)

        stop = stopping_times(t, scheme, C, m)
        N_c[block], censored[block] = censored_mle(t, stop)
        uncensored[block] = n / t.sum(axis=1)

        if C is not None:
            by_time[block] = stop == C
            n_without_time_censoring += np.count_nonzero(t.max(axis=1) < C)

    if n_without_time_censoring:
        warnings.warn(
            f"No time censoring occurred in {n_without_time_censoring} of "
            f"{n_simulations} simulations. Parameter C might be too high. C = {C}"
        )

    return SimulationResult(scheme, lam, n, C, m, uncensored, censored, N_c, by_time)


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------

def plot_standardized_mle(mles, lam, n_s, title, xlabel):
    """Histogram of sqrt(D) * (mle - lam) / mle with N(0,1) and kernel density."""
    standardized = np.sqrt(n_s) * (mles - lam) / mles
    standardized = standardized[np.isfinite(standardized)]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(standardized, bins=50, density=True, color="lightgray", edgecolor="white")

    x = np.linspace(*ax.get_xlim(), 500)
    ax.plot(x, stats.norm.pdf(x), color="red", linewidth=2, label="N(0,1)")
    ax.plot(x, stats.gaussian_kde(standardized)(x), color="blue", linewidth=2,
            label="Kernel density")

    # small p-value -> hint for non N(0, 1) distribution
    p_value = stats.kstest(standardized, "norm").pvalue
    ax.plot([], [], " ", label=f"KS p-value = {p_value:.3g}")

    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    return fig


def plot_result(res: SimulationResult):
    """Standardized uncensored and censored MLE for one simulation result."""
    fig_u = plot_standardized_mle(
        res.uncensored_mles, res.lam, res.uncensored_N,
        title=f"Standardized uncensored MLE vs. N(0,1) using t ~ Exp({res.lam})",
        xlabel=r"$\sqrt{n}\,(\hat{\lambda} - \lambda)\,/\,\hat{\lambda}$",
    )
    fig_c = plot_standardized_mle(
        res.censored_mles, res.lam, res.censored_N,
        title=(f"Standardized censored MLE ({res.label}) vs. N(0,1) "
               f"using t ~ Exp({res.lam}), {res.parameter_text}"),
        xlabel=r"$\sqrt{N}\,(\hat{\lambda} - \lambda)\,/\,\hat{\lambda}$",
    )
    return fig_u, fig_c


def plot_n_histogram(res: SimulationResult):
    """Distribution of the number of observed failures N (see module docstring)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(res.censored_N, bins=30, color="lightgray", edgecolor="white")
    ax.set_title(f"Number of observed failures N ({res.label}, {res.parameter_text})",
                 fontsize=10)
    ax.set_xlabel("N")
    ax.set_ylabel("Frequency")
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

LAMBDA = 10
N_OBS = 1000
SEED = None   # set e.g. SEED = 2026 for reproducible results

RUNS = [
    dict(scheme="type1",   n_simulations=100_000, C=0.1),
    dict(scheme="type2",   n_simulations=10_000,  m=650),
    dict(scheme="hybrid1", n_simulations=10_000,  C=0.1, m=650),
    dict(scheme="hybrid2", n_simulations=10_000,  C=0.1, m=650),
]


def main():
    rng = np.random.default_rng(SEED)

    for run in RUNS:
        res = simulate(lam=LAMBDA, n=N_OBS, rng=rng, **run)
        plot_result(res)
        if run["scheme"] in ("hybrid1", "hybrid2"):
            plot_n_histogram(res)   # cut off at m, otherwise ~ normal around ~632

    plt.show()


if __name__ == "__main__":
    main()