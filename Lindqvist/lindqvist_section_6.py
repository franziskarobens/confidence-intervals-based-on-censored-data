"""
Fiducial inference simulation for n = 2 (Section 6 example, Lindqvist)

Setting: c = 1, one failure at Y_1 = 0.5, one censoring at Y_2 = 1
=> R = 1, S = 3/2, T = t = 2/3

For a simulated pair u = (u1, u2) with order statistics u(1) <= u(2):

  bar theta(u, t) = 2*u(1)   if u(2) >= 2*u(1)
  bar theta(u, t) = u(2)     if u(2) <  2*u(1)

Fiducial distribution (5):
  theta(u, t) = 2*u(1)   conditional on  u(2) >= 2*u(1)
"""

import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from tqdm import tqdm


rng = np.random.default_rng(2026)

c_val = 1.0
t_val = 2 / 3
n_sim = 100_000

# --- Simulate 100,000 pairs of iid standard exponential variables ---------
U1 = rng.exponential(scale=1.0, size=n_sim)   # rexp(rate = 1) -> scale = 1/rate
U2 = rng.exponential(scale=1.0, size=n_sim)

u_min = np.minimum(U1, U2)
u_max = np.maximum(U1, U2)

# --- theta(u, t): unconditional piecewise solution -------------------------
# theta = 2*u(1) if u(2) >= 2*u(1), else theta = u(2)
theta_full = np.where(u_max >= 2 * u_min, 2 * u_min, u_max)

# --- Fiducial version from Eq. (5): only keep draws where u(2) >= 2*u(1) ---
fiducial_condition = u_max >= 2 * u_min
theta_fiducial = 2 * u_min[fiducial_condition]

print("Fraction of draws satisfying u(2) >= 2*u(1):",
      round(fiducial_condition.mean(), 4), "\n") # indeed P(u(2) >= 2*u(1)) = 2/3

# --- Confidence intervals: delete 2.5% at each end of the histogram -------
# np.quantile uses linear interpolation = R's default (type 7)
CI_full = np.quantile(theta_full, [0.025, 0.975])
CI_fiducial = np.quantile(theta_fiducial, [0.025, 0.975])

print(f"95% CI (theta(u, 2/3), unconditional):  "
      f"[{CI_full[0]:.3f}, {CI_full[1]:.3f}]")
print(f"95% CI (fiducial version, conditional): "
      f"[{CI_fiducial[0]:.3f}, {CI_fiducial[1]:.3f}]\n")
# -> same values as in paper

# --------------------------------------------------------------------------
# Plot: two histograms, matching Figure 2 in Lindqvist
# --------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
x_grid = np.linspace(0, 5, 500)

# Left panel: unconditional
ax = axes[0]
ax.hist(theta_full, bins=600, density=True, color="lightgray")
for v in CI_full:
    ax.axvline(v, color="red", linestyle="--", linewidth=1)

# Given: F_theta(x) = 1 - 2 *exp(- 3/2 *x) + exp(-2*x)
# Density = d/dx F(x) = 3 * exp(- 3/2 * x) - 2 * exp(-2 * x)
density_theta = 3 * np.exp(- 1.5 * x_grid) - 2 * np.exp(-2 * x_grid)

ax.plot(x_grid, density_theta, color="red", linewidth=1,
        label="Theoretical density")
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.5)
ax.set_xlabel(r"$\theta$")
ax.set_ylabel("Density")
ax.legend(loc="upper right", frameon=False)

# Right panel: fiducial (conditional)
ax = axes[1]
ax.hist(theta_fiducial, bins=600, density=True, color="lightgray")
for v in CI_fiducial:
    ax.axvline(v, color="red", linestyle="--", linewidth=2)

density_theta = 3 / 2 * np.exp(-3 / 2 * x_grid)
ax.plot(x_grid, density_theta, color="red", linewidth=1,
        label="Theoretical density Exp(3/2)")
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.5)
ax.set_xlabel(r"$\theta$")
ax.set_ylabel("Density")
ax.legend(loc="upper right", frameon=False)

plt.tight_layout()
plt.show()


# ----- Helper function to get theta using eq (4) --------------------------------

@njit
def theta_from_u(u, t, c):
    n = len(u)
    u = np.sort(u)
    # A[k] = u(1) + ... + u(k), A[0] = 0   (R: A[k + 1])
    A = np.concatenate((np.zeros(1), np.cumsum(u)))
    theta_best = 0.0

    for k in range(0, n + 1):
        # R: u[k] -> Python u[k-1];   R: u[k+1] -> Python u[k]
        if k == 0:
            lower = 0.0
        else:
            lower = u[k - 1] / c

        if k == n:
            upper = np.inf
        else:
            upper = u[k] / c

        if k == 0:
            theta_best = max(theta_best, upper)
        else:
            A_k = A[k]

            if k == n:
                candidate = t * A_k / n
                if candidate >= lower:
                    theta_best = max(theta_best, candidate)
            else:
                denominator = k - t * (n - k) * c

                tau_upper = k * upper / (A_k + (n - k) * c * upper)

                if t >= tau_upper:
                    theta_best = max(theta_best, upper)
                else:
                    if denominator > 0:
                        candidate = t * A_k / denominator
                        if lower <= candidate <= upper:
                            theta_best = max(theta_best, candidate)

    return theta_best


# ----- Check coverage probabilities of bounds and CI -----------------------------------------
# see end of page 5 in Lindqvist

n = 10
theta_true = 1.0
n_simulation = 10_000
m = 1000
alpha = 0.05

R_values = np.zeros(n_simulation)
T_values = np.zeros(n_simulation)
lower_bounds = np.zeros(n_simulation)
upper_bounds = np.zeros(n_simulation)

for sim in tqdm(range(n_simulation)):

    # Generate n = 10 exponential lifetimes with true hazard theta = 1
    # (rate = theta -> scale = 1/theta)
    X = rng.exponential(scale=1 / theta_true, size=n)

    # Type-I censoring at c
    Y = np.minimum(X, c_val)

    # Number of observed failures
    R = np.sum(X <= c_val)
    R_values[sim] = R

    # Total observed time
    TTT = Y.sum()

    # MLE: T = R / TTT
    t_hat = R / TTT
    T_values[sim] = t_hat

    # Generate m = 1000 vectors U, each row U = (U1,...,U10) with Ui ~ Exp(1)
    U_matrix = rng.exponential(scale=1.0, size=(m, n))

    # Calculate theta(U, t) for each of the m simulated U-vectors
    theta_values = np.empty(m)
    for j in range(m):
        theta_values[j] = theta_from_u(U_matrix[j, :], t_hat, c_val)

    # Fiducial/confidence limits: remove 2.5% from each tail
    bounds = np.quantile(theta_values, [alpha / 2, 1 - alpha / 2])

    lower_bounds[sim] = bounds[0]
    upper_bounds[sim] = bounds[1]


# Coverage calculations
lower_coverage = np.mean(lower_bounds <= theta_true)
upper_coverage = np.mean(upper_bounds >= theta_true)
two_sided_coverage = np.mean(
    (lower_bounds <= theta_true) & (upper_bounds >= theta_true)
)

# Results
print()
print("=" * 60)
print("Simulation results")
print("=" * 60)

print("n =", n)
print("c =", c_val)
print("true theta =", theta_true)
print("number of datasets =", n_simulation)
print("m =", m, "\n")

print("Coverage of lower bound:")
print(f"{lower_coverage}\n")

print("Coverage of upper bound:")
print(f"{upper_coverage}\n")

print("Coverage of two-sided interval:")
print(f"{two_sided_coverage}\n")