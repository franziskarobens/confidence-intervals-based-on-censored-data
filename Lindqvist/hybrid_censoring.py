"""
Plot of the data generating function (DGF) of the MLE, eq. (8),
for Hybrid Type I and Hybrid Type II censoring.

Setting: n iid lifetimes X_i ~ Exp(theta), censoring time c, r = number of
failures that stops the experiment early.

    Hybrid Type I : experiment stops at min(X_(r), c)
    Hybrid Type II: experiment stops at max(X_(r), c)

With stopping time `stop`:

    R   = number of failures up to `stop`
    S   = sum_i min(X_i, stop)        (total time on test)
    MLE = T = R / S

Refactoring idea
----------------
Eq. (8) is just the MLE evaluated at X = U / theta, where U ~ Exp(1):

    tau(U, theta) = T(U / theta)

So one function `hybrid_mle` is used both for the observed MLE T_obs (from a
sample X) and for the curve tau(U, theta). The four case distinctions of the R
script (R < r / R >= r, for both schemes) are covered by the stopping time.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

SCHEMES = {
    "hybrid1": "Hybrid Type I",
    "hybrid2": "Hybrid Type II",
}


def hybrid_mle(x, r, c, scheme):
    """
    MLE T = R / S of the exponential rate under hybrid censoring.

    `x` has shape (..., n); the MLE is computed along the last axis,
    so a whole grid of samples can be processed at once.
    """
    if scheme not in SCHEMES:
        raise ValueError(f"scheme must be one of {list(SCHEMES)}, got {scheme!r}")

    x = np.atleast_2d(x)
    n = x.shape[-1]
    if not 1 <= r <= n:
        raise ValueError(f"r must satisfy 1 <= r <= n. r = {r} and n = {n}")

    x_r = np.partition(x, r - 1, axis=-1)[..., r - 1]          # r-th order statistic
    stop = np.minimum(x_r, c) if scheme == "hybrid1" else np.maximum(x_r, c)
    stop = stop[..., None]

    R = (x <= stop).sum(axis=-1)                # observed failures
    S = np.minimum(x, stop).sum(axis=-1)        # total time on test
    return R / S


def observed_mle(n, true_theta, r, c, scheme, rng):
    """Draw one sample X ~ Exp(true_theta) and return its MLE T_obs."""
    X = rng.exponential(scale=1 / true_theta, size=n)    # rexp(rate = true_theta)
    return hybrid_mle(X, r, c, scheme)[0]


def tau_curve(U, r, c, scheme, n_points=1000):
    """tau(U, theta) = T(U / theta) on a grid theta in [0, max(U / c)]."""
    theta = np.linspace(0, np.max(U / c), n_points)
    tau = np.zeros(n_points)                  # tau(U, 0) = 0 (limit theta -> 0)
    tau[1:] = hybrid_mle(U[None, :] / theta[1:, None], r, c, scheme)
    return theta, tau


def plot_dgf_mle(n, true_theta, r, c, scheme, rng):
    """Simulate X (for T_obs) and U (for the curve) and plot tau(U, theta)."""
    T_obs = observed_mle(n, true_theta, r, c, scheme, rng)
    U = rng.exponential(scale=1.0, size=n)               # rexp(rate = 1)
    theta, tau = tau_curve(U, r, c, scheme)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(theta, tau, color="black", linewidth=1.5,
            label=r"$\tau(\mathbf{U}, \theta)$")
    for i, u in enumerate(U / c):
        ax.axvline(u, color="gray", linestyle=":", linewidth=1,
                   label=r"$U_i / c$" if i == 0 else None)
    ax.axhline(T_obs, color="red", label=r"$T_{obs}$")

    ax.set_title(f"DGF of the MLE, {SCHEMES[scheme]} censoring "
                 f"(n = {n}, r = {r}, c = {c})", fontsize=10)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\tau$")
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    return fig, T_obs, U


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

SEED = None   # set e.g. SEED = 2026 for reproducible results

RUNS = [
    dict(scheme="hybrid2", n=8,  true_theta=1, c=1, r=3),
    dict(scheme="hybrid1", n=10, true_theta=1, c=1, r=5),
]


def main():
    rng = np.random.default_rng(SEED)
    for run in RUNS:
        plot_dgf_mle(rng=rng, **run)
    plt.show()


if __name__ == "__main__":
    main()