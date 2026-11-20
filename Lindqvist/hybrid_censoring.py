"""
Plot of the DGF of the MLE, eq. (8),
for Hybrid Type I and Hybrid Type II censoring.

Setting: n iid lifetimes X_i ~ Exp(theta), censoring time c, r = number of
failures that stops the experiment early.

    Hybrid Type I : experiment stops at min(X_(r), c)
    Hybrid Type II: experiment stops at max(X_(r), c)

With stopping time `stop`:

    R   = number of failures up to `stop`
    TTT   = sum_i min(X_i, stop)        (total time on test)
    MLE = T = R / TTT
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def hybrid_mle(x, r, c, scheme):
    """
    MLE: T = R / TTT of the exponential rate under hybrid censoring.
    """
    if scheme not in {"hybrid1", "hybrid2"}:
        raise ValueError(f"scheme must be one of {list({"hybrid1", "hybrid2"})}, got {scheme!r}")

    n = len(x)
    if not 1 <= r <= n:
        raise ValueError(f"r must satisfy 1 <= r <= n. r = {r} and n = {n}")

    x_r = np.sort(x)[r-1]                         # r-th order statistic
    stop = np.minimum(x_r, c) if scheme == "hybrid1" else np.maximum(x_r, c)

    R = sum(x <= stop)                            # observed failures
    TTT = np.minimum(x, stop).sum()               # total time on test
    return R / TTT


def tau_curve(U, r, c, scheme, n_points=1000):
    """tau(U, theta) = T(U / theta) on a grid theta in [0, max(U / c)]."""
    theta = np.linspace(start=0, stop=np.max(U / c), num=n_points)
    tau = [0]                  # tau(U, 0) = 0 (limit theta -> 0)
    for th in theta[1:]:
        tau.append(hybrid_mle(U / th, r, c, scheme))
    return theta, tau


def plot_dgf_mle(n, true_theta, r, c, scheme, rng):
    """Simulate X (for T_obs) and U (for the curve) and plot tau(U, theta)."""

    U = rng.exponential(scale=1.0, size=n)               # rexp(rate = 1)
    theta, tau = tau_curve(U, r, c, scheme)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(theta, tau, color="black", linewidth=1.5,
            label=r"$\tau(\mathbf{U}, \theta)$")
    for i, u in enumerate(U / c):
        ax.axvline(u, color="gray", linestyle=":", linewidth=1,
                   label=r"$U_i / c$" if i == 0 else None)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\tau$")
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    return fig, U


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    plot_dgf_mle(n=6, true_theta=1, r=3, c=1, scheme="hybrid2", rng=rng)