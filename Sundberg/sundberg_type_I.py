"""For type I censored data the fisher info is given by
I(theta) = n * (1 - exp(- theta * c)) / theta^2
which delivers that 
1 - alpha = Pr(-z <= sqrt(n * (1 - exp(- theta * c))) * (1 - theta / hat theta)^2 <= z) 
Solve this inequality numerically and construct a CI.
"""

import numpy as np
from scipy.stats import norm
from collections import namedtuple
from scipy.optimize import brentq, minimize_scalar
from scipy.special import lambertw
import sundberg


def get_MLE(X, c, theta_is_rate=True):
    TTT = np.minimum(X, c).sum()
    R = np.sum(X < c)
    MLE = R / TTT if theta_is_rate else TTT / R
    return TTT, R, MLE


# Since these intervals don't have a form [lower * MLE, upper * MLE],
# we return the bounds itself, not the constants (as we did in sundberg.py)
Limits = namedtuple("Limits", ["lower", "upper"])


def C1(n, c, X, alpha):
    _, _, MLE = get_MLE(X, c)
    z = norm.ppf(1 - alpha / 2)
    term = z / np.sqrt(n * (1 - np.exp(- c * MLE))) # TODO: man hat zwei Mglk: schätzer benutzen r = n * (1 - np.exp(- c * MLE)) oder die beobachtete Zahl R einsetzen. Vergleiche beides
    return Limits(MLE * (1 - term), MLE * (1 + term))


# Helper function for C2()
def Q2(theta, MLE, n, c):
    """Pivot Q(theta) = sqrt(n (1 - exp(-c theta))) * (MLE/theta - 1)."""
    return np.sqrt(n * (-np.expm1(-c * theta))) * (MLE / theta - 1.0)


def C2(n, c, X, alpha):
    z = norm.ppf(1 - alpha / 2)   # z_{1 - alpha/2}
    _, _, MLE = get_MLE(X, c)

    f = lambda th, z: Q2(th, MLE, n, c) - z

    # lower bound: Q(lower) = z, lower in (0, MLE)
    eps = MLE * 1e-12
    lower = brentq(f, eps, MLE, args=(z,), xtol=1e-14, rtol=1e-12)

    # upper bound: Q(upper) = -z, upper in (MLE, inf)
    if -z <= -np.sqrt(n):
        upper = np.inf  # Q fällt nie unter -sqrt(n) -> keine Nullstelle
    else:
        hi = 2 * MLE
        while Q2(hi, MLE, n, c) > -z:
            hi *= 2
            if hi > 1e300:
                hi = np.inf
                break
        upper = (np.inf if np.isinf(hi)
                   else brentq(f, MLE, hi, args=(-z,), xtol=1e-14, rtol=1e-12))

    return Limits(lower, upper)


# Helper function for C3()
def Q3(theta, MLE, n, c):
    """Pivot Q3(theta) = sqrt(n (1 - exp(-c theta))) * log(MLE/theta)"""
    return np.sqrt(n * (-np.expm1(-c * theta))) * np.log(MLE / theta)


def C3(n, c, X, alpha):
    z = norm.ppf(1 - alpha / 2)
    _, _, MLE = get_MLE(X, c)

    f = lambda th, z: Q3(th, MLE, n, c) - z

    # Find maximum of Q3 on (0, MLE)
    maximum = minimize_scalar(
        lambda th: -Q3(th, MLE, n, c),
        bounds=(MLE * 1e-10, MLE),
        method="bounded"
    )

    theta_max = maximum.x
    q_max = Q3(theta_max, MLE, n, c)

    # Lower bound: Q3(lower) = z
    if q_max <= z:
        lower = 0.0
    else:
        lower = brentq(
            f,
            theta_max,
            MLE,
            args=(z,),
            xtol=1e-14,
            rtol=1e-12
        )

    # Upper bound: Q3(upper) = -z
    hi = 2 * MLE

    while Q3(hi, MLE, n, c) > -z:
        hi *= 2
        if hi > 1e300:
            hi = np.inf
            break

    upper = (
        np.inf if np.isinf(hi)
        else brentq(
            f,
            MLE,
            hi,
            args=(-z,),
            xtol=1e-14,
            rtol=1e-12
        )
    )

    return Limits(lower, upper)



# Helper function for C4()
def Q4(theta, MLE, n, c):
    """Pivot Q4(theta) = 3 * sqrt(n * (1 - exp(-c theta))) * (theta^(1/3) - MLE^(1/3)) / theta^(1/3)"""
    return (
        3 * np.sqrt(n * (-np.expm1(-c * theta)))
        * (1 - (MLE / theta)**(1/3))
    )


def C4(n, c, X, alpha):
    z = norm.ppf(1 - alpha / 2)
    _, _, MLE = get_MLE(X, c)

    # Special case: MLE = 0
    if MLE == 0:
        lower = 0.0

        if z >= 3 * np.sqrt(n):
            upper = np.inf
        else:
            f = lambda th: Q4(th, MLE, n, c) - z

            hi = 2.0 / c
            while f(hi) < 0:
                hi *= 2

            upper = brentq(
                f, 0.0 + 1e-14, hi,
                xtol=1e-14, rtol=1e-12
            )

        return Limits(lower, upper)

    # Find minimum of Q4 below MLE
    minimum = minimize_scalar(
        lambda th: Q4(th, MLE, n, c),
        bounds=(MLE * 1e-10, MLE),
        method="bounded"
    )

    theta_min = minimum.x
    q_min = Q4(theta_min, MLE, n, c)

    # Lower bound: Q4(lower) = -z
    if q_min >= -z:
        lower = 0.0
    else:
        lower = brentq(
            lambda th: Q4(th, MLE, n, c) + z,
            theta_min,
            MLE,
            xtol=1e-14,
            rtol=1e-12
        )

    # Upper bound: Q4(upper) = z
    if z >= 3 * np.sqrt(n):
        upper = np.inf
    else:
        hi = 2 * MLE

        while Q4(hi, MLE, n, c) < z:
            hi *= 2

        upper = brentq(
            lambda th: Q4(th, MLE, n, c) - z,
            MLE,
            hi,
            xtol=1e-14,
            rtol=1e-12
        )

    return Limits(lower, upper)


def C5(n, c, X, alpha):
    _, R, MLE = get_MLE(X, c)
    if R == 0:
        return None                       # MLE not defined
    constants = sundberg.C5(r=R, alpha=alpha)
    return Limits(constants.c_lower * MLE, constants.c_upper * MLE)


# check the coverage probability with monte carlo simulation
if __name__ == "__main__":

    n_sim = 10_000
    alpha = 0.01
    theta = 10
    c = 0.5 * theta
    n = 100

    coverage = 0

    for i in range(n_sim):
        # Create data
        rng = np.random.default_rng(i)
        X = rng.exponential(scale=1/theta, size=n)
     
        lower, upper = C5(n, c, X, alpha)

        if lower <= theta <= upper:
            coverage += 1

    print(1 - coverage / n_sim) # should be alpha

