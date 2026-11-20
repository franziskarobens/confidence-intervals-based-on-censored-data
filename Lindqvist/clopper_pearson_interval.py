import numpy as np
from scipy.stats import beta

# --- Example: Clopper-Pearson Interval (Section 2.1 in Lindqvist) --------------


def clopper_pearson_interval(T, n, alpha, kind="two-sided"):
    """kind: "upper", "lower" or "two-sided". """
    if kind == "upper":
        lower_bound = 0.0
        upper_bound = 1.0 if T == n else beta.ppf(1 - alpha, T + 1, n - T)
    elif kind == "lower":
        lower_bound = 0.0 if T == 0 else beta.ppf(alpha, T, n + 1 - T)
        upper_bound = 1.0
    elif kind == "two-sided":
        tail_alpha = alpha / 2
        lower_bound = 0.0 if T == 0 else beta.ppf(tail_alpha, T, n + 1 - T)
        upper_bound = 1.0 if T == n else beta.ppf(1 - tail_alpha, T + 1, n - T)

    return lower_bound, upper_bound


if __name__ == "__main__":
    theta = 0.5 # in [0, 1]
    n = 100
    alpha = 0.01
    n_sim = 100_00
    rng = np.random.default_rng()
    T = rng.binomial(n=n, p=theta, size=n_sim)

    coverage = 0
    for T_obs in T:
        lower_bound, upper_bound = clopper_pearson_interval(T_obs, n, alpha=alpha, kind="two-sided")
        print(lower_bound, upper_bound)
        if lower_bound <= theta <= upper_bound:
            coverage += 1
    print(coverage / n_sim)




