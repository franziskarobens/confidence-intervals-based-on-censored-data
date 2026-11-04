import warnings
from collections import namedtuple

import numpy as np
from scipy.special import lambertw
from scipy.stats import chi2, norm

# constants c_lower and c_upper from C1, ..., C7 from Sundberg
# assuming parametrization is rate = theta

Limits = namedtuple("Limits", ["c_lower", "c_upper"])


def C1(r, alpha):
    z = norm.ppf(1 - alpha / 2)
    r = np.asarray(r, dtype=float)
    return Limits(1 - z / np.sqrt(r),
                  1 + z / np.sqrt(r))


def C2(r, alpha):
    z = norm.ppf(1 - alpha / 2)
    r = np.asarray(r, dtype=float)
    if np.any(r < z**2):
        warnings.warn("C2 limits aren't valid. need higher r.")
    return Limits(1 / (1 + z / np.sqrt(r)),
                  1 / (1 - z / np.sqrt(r)))


def C3(r, alpha):
    z = norm.ppf(1 - alpha / 2)
    r = np.asarray(r, dtype=float)
    return Limits(np.exp(-z / np.sqrt(r)),
                  np.exp(z / np.sqrt(r)))


def C4(r, alpha):
    z = norm.ppf(1 - alpha / 2)
    r = np.asarray(r, dtype=float)
    return Limits(1 / (1 + z / (3 * np.sqrt(r)))**3,
                  1 / (1 - z / (3 * np.sqrt(r)))**3)


def C5(r, alpha):
    z = norm.ppf(1 - alpha / 2)
    r = np.asarray(r, dtype=float)
    t = z**2 / (2 * r)

    arg = -np.exp(-t - 1)              # lies in (-1/e, 0) for t > 0
    c1 = -lambertw(arg, k=0).real      # solution in (0, 1)
    c2 = -lambertw(arg, k=-1).real     # solution in (1, inf)
    return Limits(c1, c2)


def C6(r, alpha):
    r = np.asarray(r, dtype=float)
    c1 = chi2.ppf(alpha / 2, df=2 * r)
    c2 = chi2.ppf(1 - alpha / 2, df=2 * r)
    return Limits(c1 / (2 * r), c2 / (2 * r))


def C7(r, alpha):
    r = np.asarray(r, dtype=float)
    c1 = chi2.ppf(alpha / 2, df=2 * r + 1)
    c2 = chi2.ppf(1 - alpha / 2, df=2 * r + 1)
    return Limits(c1 / (2 * r), c2 / (2 * r))


# combine all Ci methods from Sundberg in one function
def CI(method_nr, r, alpha, theta_is_rate=True):
    match method_nr:
        case 1: c1, c2 = C1(r, alpha)
        case 2: c1, c2 = C2(r, alpha)
        case 3: c1, c2 = C3(r, alpha)
        case 4: c1, c2 = C4(r, alpha)
        case 5: c1, c2 = C5(r, alpha)
        case 6: c1, c2 = C6(r, alpha)
        case 7: c1, c2 = C7(r, alpha)
        case _: raise ValueError(f"Unknown Method: {method_nr}, method number must be in 1, ..., 7")

    if theta_is_rate:
        return c1, c2
    return 1 / c2, 1 / c1