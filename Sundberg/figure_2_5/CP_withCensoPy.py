"""Compute non-coverage probabilities for all combinations of the parameter lists
below and save the results in "non_coverage_probs".
nu = np.inf means no censoring.
ci_type = 1: two CSVs per combination (lower-bounded CI, upper-bounded CI)
ci_type = 2: one CSV per combination (two-sided CI)"""
import sys
from pathlib import Path

# parent folder of the current working directory (notebook)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> ...\Sundberg

import numpy as np
import pandas as pd
from CensoPy import *
import sundberg
from tqdm import tqdm
from scipy.stats import chi2
from itertools import product

# ---------------------------------------------------------------------------
# Parameters: every list is run through, all combinations are computed.
# ---------------------------------------------------------------------------
theta_list = [1]
nu_list = [0.5]    # in {0.5, 1, 2, np.inf}; np.inf = no censoring
alpha_list = [0.1]              # in {0.01, 0.05, 0.1}
Digits_list = [100]              # works for >= 50 (ignored if nu = np.inf)
ci_type = 2                      # 1: one sided CIs, 2: two sided CI

n_range = range(4, 51)           # 4, ..., 50
methods = range(1, 8)            # C1, ..., C7

overwrite = True                 # False: skip combinations whose CSV already exists

out_dir = Path(__file__).resolve().parent
csv_dir = out_dir / "non_coverage_probs"
csv_dir.mkdir(exist_ok=True)


def make_cdf_mle(n, theta, nu, Digits):
    # Return s -> P(MLE <= theta / s) for the given censoring scheme.
    if np.isinf(nu):
        # no censoring: 2 * sum(X) / theta ~ chi^2_{2n}, MLE = sum(X) / n
        return lambda s: chi2.cdf(2 * n / s, df=2 * n)

    # does not depend on the method, so create it once per n
    c = nu * theta
    HCS = C_I_HCS(n=n,
                  thresholds=[c],
                  censoring_numbers=[n],
                  Digits=Digits,
                  Method="analytic")
    return lambda s: HCS.cdf_MLE(t=theta / s, theta=theta)


def kinds_for(ci_type):
    if ci_type == 1:
        return ["lower", "upper"]
    if ci_type == 2:
        return ["two_sided"]
    raise ValueError(f"Unknown ci_type: {ci_type}")


def compute_non_coverage(theta, nu, alpha, Digits, ci_type):
    """Return a dict {kind: DataFrame}."""
    kinds = kinds_for(ci_type)
    non_coverage = {kind: {m: [] for m in methods} for kind in kinds}

    # one-sided CI of level alpha = two-sided CI of level 2 * alpha, one tail each
    alpha_ci = alpha if ci_type == 2 else 2 * alpha

    for n in tqdm(n_range):
        # R ~ Bin(n, 1 - exp(-nu)), use its expectation for r
        # (for nu = inf this gives r = n)
        r = n * (1 - np.exp(-nu))

        cdf_mle = make_cdf_mle(n, theta, nu, Digits)

        for m in methods:
            # use rate = 1 / theta bc Censopy uses this parametrisation
            c_lower, c_upper = sundberg.CI(method_nr=m, r=r, alpha=alpha_ci,
                                           theta_is_rate=False)

            # C2 with r < z^2: c_lower <= 0, so the event cannot happen
            if c_lower > 0 and np.isfinite(c_lower):
                prob_lower = cdf_mle(c_lower)
            else:
                prob_lower = 1.0

            if c_upper > 0 and np.isfinite(c_upper):
                prob_upper = cdf_mle(c_upper)
            else:
                prob_upper = 0.0

            if ci_type == 1:
                non_coverage["lower"][m].append(1 - prob_lower)  # NCP of lower bounded CI
                non_coverage["upper"][m].append(prob_upper)      # NCP of upper bounded CI
            else:
                non_coverage["two_sided"][m].append(1 - prob_lower + prob_upper)

    dfs = {}
    for kind in kinds:
        df = pd.DataFrame({"n": list(n_range)})
        for m in methods:
            df[f"C{m}"] = non_coverage[kind][m]
        dfs[kind] = df
    return dfs


combinations = list(product(theta_list, nu_list, alpha_list, Digits_list))

for i, (theta, nu, alpha, Digits) in enumerate(combinations, start=1):
    label = f"[{i}/{len(combinations)}] theta={theta}, nu={nu}, alpha={alpha}, Digits={Digits}"

    csv_files = {
        kind: csv_dir / f"non_coverage_{kind}_nu{nu}_alpha{alpha}_theta{theta}_Digits{Digits}.csv"
        for kind in kinds_for(ci_type)
    }

    if not overwrite and all(f.exists() for f in csv_files.values()):
        print(f"{label}: skipped, files exist")
        continue

    print(f"{label}: running")
    dfs = compute_non_coverage(theta, nu, alpha, Digits, ci_type)
    for kind, df in dfs.items():
        df.to_csv(csv_files[kind], index=False)   # saved right after each combination
        print(f"{label}: saved {csv_files[kind].name}")