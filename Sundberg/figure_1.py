# Replicate Figure 1a and 1b in Sundberg

import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sundberg


alpha = 0.025
N_range = range(1, 26)
CI_nrs = range(1, 8)

rows = []
for N in N_range:
    for nr in CI_nrs:
        try:
            with warnings.catch_warnings(), np.errstate(divide="ignore", invalid="ignore"):
                warnings.simplefilter("ignore")   # C2 warns for r < z^2
                c_lower, c_upper = sundberg.CI(method_nr=nr, r=N, alpha=alpha,
                                               theta_is_rate=False) # set 'theta_is_rate=False' to get exact results from Sundberg
            c_lower, c_upper = float(c_lower), float(c_upper)
        except Exception as e:
            warnings.warn(f"C{nr}, N = {N}: {e}")
            c_lower, c_upper = np.nan, np.nan
        rows.append((N, nr, c_lower, c_upper))

res = pd.DataFrame(rows, columns=["N", "CI_nr", "c_lower", "c_upper"])

linestyles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (1, 1))]


def plot_c(col, ylim, ylabel, filename):
    fig, ax = plt.subplots(figsize=(6, 4.5))

    for nr in CI_nrs:
        sub = res[res["CI_nr"] == nr]
        ax.plot(sub["N"], sub[col], color=f"C{nr - 1}",
                linestyle=linestyles[nr - 1], linewidth=1.5, label=f"C{nr}")

    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(ylabel)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.3)

    ax.legend(title="Method", loc="upper center", bbox_to_anchor=(0.5, -0.15),
              ncol=7, frameon=False, columnspacing=1.0, handlelength=2.5)

    fig.tight_layout()
    fig.savefig(filename, bbox_inches="tight")
    return fig


plot_c("c_lower", (0, 0.8), r"$c_l$", "fig_1a_sundberg_c_lower.pdf")
plot_c("c_upper", (1, 8), r"$c_u$", "fig_1b_sundberg_c_upper.pdf")
plt.show()

# C1 and C2 differ from the paper because we use another parametrisation
# (rate = theta instead of rate = 1 / theta)