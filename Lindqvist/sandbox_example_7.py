"""Implement my example:

Suppose T ∼ Exp(1 / θ) and U ∼ Exp (1).
τ(U, θ) = θ · U
"""

import numpy as np
import matplotlib.pyplot as plt


def a(t, alpha):
    return -t / np.log(1 - alpha)


def sim_theorem_1(alpha=0.1, theta=1.0, n_sim=100_000, n_min=4, n_max=50, seed=None):
    """Simulate coverage of the upper CB a(T) for T ~ Exp(1/theta), i.e. mean theta."""
    rng = np.random.default_rng(seed)
    n_range = np.arange(n_min, n_max + 1)

    coverages = []
    for _ in n_range:
        T = rng.exponential(scale=theta, size=n_sim)  # scale = mean = theta
        coverages.append(np.mean(a(T, alpha) >= theta))
    coverages = np.array(coverages)

    _, ax = plt.subplots()
    ax.plot(n_range, coverages, "o-", label="simulated coverage")
    ax.axhline(1 - alpha, color="red", linestyle="--", label=f"nominal 1 - α = {1 - alpha}")
    ax.set_xlabel("n")
    ax.set_ylabel("coverage of a(T)")
    ax.set_ylim([0.8, 1])
    ax.set_title(f"Coverage of upper CB (α = {alpha}, {n_sim:,} simulations)")
    ax.legend()
    plt.show()

    return n_range, coverages


if __name__ == "__main__":
    sim_theorem_1()

