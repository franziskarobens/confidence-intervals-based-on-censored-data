"""Read saved non-coverage probabilities (CSV) and plot them.

Put this file in the same folder as the script that creates the data
(the folder containing "non_coverage_probs").

ci_type = 1: two figures (lower-bounded CI, upper-bounded CI)
ci_type = 2: one figure (two-sided CI)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Parameters: must match the values used when the data was created, because
# they are part of the file name. Write them exactly as in the original script
# (e.g. theta = 1 gives "theta1", theta = 1.0 gives "theta1.0").
# ---------------------------------------------------------------------------
nu = 0.5                  # in {0.5, 1, 2, np.inf}
alpha = 0.1                 # in {0.01, 0.05, 0.1}; for ci_type = 1 the one-sided level
theta = 1
Digits = 100
ci_type = 2                  # 1: one sided CI, 2: two sided CI

methods = range(1, 8)        # C1, ..., C7  (use e.g. [5, 6, 7] for a subset)
y_limits = [0, 0.2]        # y-axis shows y_limits; None = automatic
save_plot = True             # save the figures as png and pdf

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
base_dir = Path(__file__).resolve().parent
csv_dir = base_dir / "non_coverage_probs"
png_dir = base_dir / "non_coverage_plots" / "png"
pdf_dir = base_dir / "non_coverage_plots" / "pdf"

if ci_type == 1:
    kinds = ["lower", "upper"]
elif ci_type == 2:
    kinds = ["two_sided"]
else:
    raise ValueError(f"Unknown ci_type: {ci_type}")

kind_titles = {
    "lower": "lower-bounded CI",
    "upper": "upper-bounded CI",
    "two_sided": "two-sided CI",
}


def plot_kind(kind):
    stem = f"non_coverage_{kind}_nu{nu}_alpha{alpha}_theta{theta}_Digits{Digits}"
    csv_file = csv_dir / f"{stem}.csv"

    if not csv_file.exists():
        available = "\n  ".join(p.name for p in sorted(csv_dir.glob("*.csv")))
        raise FileNotFoundError(
            f"{csv_file} not found.\nAvailable files:\n  {available}"
        )

    df = pd.read_csv(csv_file)

    fig, ax = plt.subplots(figsize=(8, 5))

    for m in methods:
        ax.plot(df["n"], df[f"C{m}"], marker="o", markersize=3, label=f"C{m}")

    ax.axhline(alpha, linestyle="--", color="grey", label=r"nominal $\alpha$")
    ax.set_xlabel("n")
    ax.set_ylabel("non-coverage probability")
    ax.set_title(rf"{kind_titles[kind]}: $\nu$ = {nu:g}, $\alpha$ = {alpha:g}, "
                 rf"$\theta$ = {theta:g}, Digits = {Digits}")

    if y_limits is not None:
        ax.set_ylim(y_limits)

    ax.legend(ncol=2)

    if save_plot:
        png_dir.mkdir(parents=True, exist_ok=True)
        pdf_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
        fig.savefig(pdf_dir / f"{stem}.pdf", bbox_inches="tight")

    return fig


for kind in kinds:
    plot_kind(kind)

plt.show()