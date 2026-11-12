"""Use Theorem 2, Corollary 3 and Corollary 4 from Lindqvist"""

from IPython.core.completer import __main__
import numpy as np

# compute observation T artificially (T = R / TTT for rate = theta)
def compute_tau(theta, U, c, theta_is_rate=True):
    rate = theta if theta_is_rate else 1 / theta
    numerator = np.sum(U / rate < c)
    denominator = np.sum(np.minimum(U / rate, c))
    return numerator / denominator

# uses parametrisation rate = theta
def solve_theta_for_u(U, t, c):
    # use (A1) in the Appendix A from Lindqvist
    # (loop index i is kept 1-based as in the paper
    #  array access therefore uses i-1)

    sorted_U = np.sort(U)
    n = len(U)

    green_bound = 1 / (n * c)
    yellow_bounds = np.zeros(n - 1)
    orange_bounds = np.zeros(n)

    for i in range(1, n + 1):
        numerator = i * sorted_U[i - 1]
        denominator = sorted_U[:i].sum() + (n - i) * sorted_U[i - 1]
        denominator = denominator * c
        orange_bounds[i - 1] = numerator / denominator

    for i in range(1, n):
        numerator = i * sorted_U[i]
        denominator = sorted_U[:i + 1].sum() + (n - i - 1) * sorted_U[i]
        denominator = denominator * c
        yellow_bounds[i - 1] = numerator / denominator

    pink_bounds = np.concatenate(([np.nan], yellow_bounds))

    # Case A
    if 0 <= t <= green_bound:
        return sorted_U[0] / c

    # Case B
    for i in range(2, n + 1):
        if pink_bounds[i - 1] < t <= orange_bounds[i - 1]:
            return sorted_U[i - 1] / c

    # Case C
    for i in range(1, n):
        if orange_bounds[i - 1] < t <= yellow_bounds[i - 1]:
            return t * sorted_U[:i].sum() / (i - (n - i) * c * t)

    # Case D
    if t > orange_bounds[n - 1]:
        return t * sorted_U.sum() / n

    raise RuntimeError("No result for solve_theta_for_u()")


# -------- Use Theorem 2/ Corollary 3 / 4 from Lindqvist --------------

# Precondition: DGF must be strictly increasing and bar theta has continuous dist'n
# This implementation only works for exponentially distributed data

# Params:
# T: observed value of t
# n: length of the vectors U_i
# c: censoring time

# k = 1,  m = 199  -> 0.99 CI (both sided)
# k = 5,  m = 999  -> 0.99 CI (both sided)
# k = 5,  m = 99   -> 0.9  CI (both sided)
# k = 50, m = 999  -> 0.9  CI (both sided)

# lower and upper are 1 - k/(m+1) CBs
# (lower, upper) is a 1 - 2k/(m+1) CI
def CI_corollary_4(T, n, c, k, m=100):

    if T == 0:
        raise ValueError("We condition on t > 0 to use Corollary 4.")

    # draw m samples of U
    bar_theta = np.zeros(m)

    for i in range(m):
        rng = np.random.default_rng(i)
        U_sim = rng.exponential(scale=1.0, size=n)
        # solve for theta
        bar_theta[i] = solve_theta_for_u(U_sim, t=T, c=c)

    ordered_theta = np.sort(bar_theta)

    lower = ordered_theta[k - 1]       # R: sorted_theta[k]
    upper = ordered_theta[m - k]       # R: sorted_theta[m - k + 1]

    #print(f"\nResults for k = {k} and m = {m}:")
    #print(f"Lower bound for theta with alpha = {round(1 - k/(m+1), 3)}   :   {round(lower, 3)}")
    #print(f"Upper bound for theta with alpha = {round(1 - k/(m+1), 3)}   :   {round(upper, 3)}")
    #print(f"CI for theta with alpha = {round(1 - 2*k/(m+1), 3)}   :   [{round(lower, 3)}, {round(upper, 3)}]")

    return lower, upper


if __name__ == "__main__":

    alpha = 0.1
    n = 100
    n_sim = 1_00
    theta = 1
    nu = 0.5
    c = nu * theta

    covered = 0

    for i in range(n_sim):
        if i % 10 == 0:
            print(f"n_sim = {i}")
        rng = np.random.default_rng(i)
        U = rng.exponential(scale=1.0, size=n)
        T = compute_tau(theta, U, c)
        CI_corollary_4(T, n, c, k=5, m = 100)

        lower, upper = CI_corollary_4(T, n, c, k=5, m = 99) # k=5, m=99 for both sided 90% interval
        covered += lower <= theta <= upper


    print(covered / n_sim)   # should be close to 1 - alpha