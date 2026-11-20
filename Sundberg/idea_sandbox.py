"""
Compare two plug-in variants of C1 for type I censored exponential data:

  (A) sundberg_type_I.C1: Varianz-Plug-in mit n * (1 - exp(-c * MLE))
  (B) C1_R:               Varianz-Plug-in mit der beobachteten Anzahl R

- plot n vs P(non coverage) for both methods (simulation)
- plot n vs mean width of the interval
"""
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

import sundberg_type_I
from sundberg_type_I import Limits


def C1_R(n, c, X, alpha):
    z = norm.ppf(1 - alpha / 2)
    R = np.sum(X < c)
    TTT = np.minimum(X, c).sum()
    MLE = R / TTT                      # Rate
    term = z / np.sqrt(R)
    return Limits(MLE * (1 - term), MLE * (1 + term))


alpha = 0.1
theta = 1.0                            # Rate
nu = 0.5
c = nu * theta

n_range = list(range(4, 51))
n_sim = 10                          # bei Bedarf auf 100_000 erhöhen (dauert deutlich länger)

rng = np.random.default_rng(2026)

noncov_A, noncov_B = [], []
width_A, width_B = [], []

for n in n_range:
    print(n, end=" ", flush=True)
    miss_A = miss_B = 0
    w_A, w_B = [], []
    n_valid = 0

    for _ in range(n_sim):
        X = rng.exponential(scale=1 / theta, size=n)
        R = np.sum(X < c)

        if R == 0:                     # MLE = 0, Intervall nicht definiert
            continue
        n_valid += 1

        TTT = np.minimum(X, c).sum()
        MLE = R / TTT                  # Rate-MLE

        C_A = sundberg_type_I.C1(n, c, MLE, alpha)
        C_B = C1_R(n, c, X, alpha)

        if not (C_A.lower <= theta <= C_A.upper):
            miss_A += 1
        if not (C_B.lower <= theta <= C_B.upper):
            miss_B += 1

        w_A.append(C_A.upper - C_A.lower)
        w_B.append(C_B.upper - C_B.lower)

    noncov_A.append(miss_A / n_valid)
    noncov_B.append(miss_B / n_valid)
    width_A.append(np.mean(w_A))
    width_B.append(np.mean(w_B))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))

ax[0].plot(n_range, noncov_A, color="red", label=r"C1 Plug-in $n(1-e^{-c\hat\theta})$")
ax[0].plot(n_range, noncov_B, color="blue", label="C1 Plug-in $R$ (beobachtet)")
ax[0].axhline(alpha, color="gray", ls="--", label=r"nominal $\alpha$")
ax[0].set_xlabel("n")
ax[0].set_ylabel("P(non coverage)")
ax[0].legend()

ax[1].plot(n_range, width_A, color="red", label=r"C1 Plug-in $n(1-e^{-c\hat\theta})$")
ax[1].plot(n_range, width_B, color="blue", label="C1 Plug-in $R$ (beobachtet)")
ax[1].set_xlabel("n")
ax[1].set_ylabel("mean width")
ax[1].legend()

plt.tight_layout()
plt.show()

print(f"{'n':>4} | {'width A':>10} | {'width B':>10} | {'Diff (A-B)':>10}")
print("-" * 44)
for n, wa, wb in zip(n_range, width_A, width_B):
    print(f"{n:>4} | {wa:>10.4f} | {wb:>10.4f} | {wa - wb:>10.4f}")


# warum habe ich immer genau dieselben weiten? 
"""


from CensoPy import *
import numpy as np
import pandas as pd
import sundberg
import matplotlib.pyplot as plt
from scipy.stats import binom


# Plot F^MLE(theta / c_l)




def plot_ruggedness():

    theta = 1
    nu = 0.5
    alpha = 0.1
    n_range = range(4, 51)
    c = nu * theta
    Digits = 50

    F_l, F_u, F_g = [], [], []

    for n in n_range:
    
        HCS = C_I_HCS(n=n,
            thresholds=[c],
            censoring_numbers=[n],
            Digits=Digits,
            Method="analytic")

        r = n * (1 - np.exp(- nu))
        c_lower, c_upper = sundberg.CI(method_nr=1, r=r, alpha=alpha, theta_is_rate=False)
        print(f"n = {n}, r = {r}:        l = {round(theta / c_lower, 2)}, u = {round(theta / c_upper, 2)}", end = " | ")
        l = - HCS.cdf_MLE(t=theta / c_lower, theta=theta)
        u = HCS.cdf_MLE(t=theta / c_upper, theta=theta)
        F_l.append(l)
        F_u.append(u)
        F_g.append(l + u)
        print(round(l + u, 2))

    non_cov_prob = 1 + np.array(F_g)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(n_range, F_l, color="red", marker="o", ms=3, label=r"$- F^{MLE}(\theta / c_{lower})$")
    ax.plot(n_range, F_u, color="blue", marker="o", ms=3, label=r"$F^{MLE}(\theta / c_{upper})$")
    ax.plot(n_range, F_g, color="green", marker="o", ms=3, label=r"$- F^{MLE}(\theta / c_{lower}) + F^{MLE}(\theta / c_{upper})$")
    ax.plot(n_range, non_cov_prob, color="black", marker="o", ms=3, label=r"$P(non coverage)$")

    # grün soll von unten gegen - alpha konvergieren
    # black soll von oben gegen alpha konvergieren

    ax.set_xlabel("r")
    ax.set_ylabel("CDFs of MLE")
    ax.legend()
    plt.tight_layout()
    plt.show()



# Plot F^MLE (theta / c_lower)
def plot_lower_MLE_CDF():
    theta = 1
    nu = 0.5
    alpha = 0.1
    n_range = range(4, 51)
    c = nu * theta
    Digits = 100

    F_l, r_range = [], []

    for n in n_range:
        print(n)
    
        HCS = C_I_HCS(n=n,
            thresholds=[c],
            censoring_numbers=[n],
            Digits=Digits,
            Method="analytic")

        r = n * (1 - np.exp(- nu))
        r_range.append(r)
        c_lower, _ = sundberg.CI(method_nr=1, r=r, alpha=alpha, theta_is_rate=False)
        F_l.append(HCS.cdf_MLE(t=theta / c_lower, theta=theta))


    fig, ax = plt.subplots(1, 3, figsize=(17, 4.5))

    ax[0].plot(n_range, F_l, color="red", marker="o", ms=3, label=r"$F^{MLE}(\theta / c_{lower})$")
    ax[0].set_xlabel("n")
    ax[0].set_ylabel("CDFs of MLE")
    ax[0].legend()

    ax[1].plot(r_range, F_l, color="red", marker="o", ms=3, label=r"$F^{MLE}(\theta / c_{lower})$")
    ax[1].set_xlabel(r"$r = n\,(1 - e^{-\nu})$")
    ax[1].set_ylabel("CDFs of MLE")
    ax[1].legend()

    ax[2].plot(n_range, r_range, color="purple", marker="o", ms=3, label=r"$r = n\,(1 - e^{-\nu})$")
    ax[2].set_xlabel("n")
    ax[2].set_ylabel("r")
    ax[2].legend()

    plt.tight_layout()
    plt.show()

# festes t, nur Einfluss von F^MLE anschauen
def plot_influence_of_CDF_MLE():
    theta = 1
    nu = 0.5
    c = nu * theta

    n_range = range(4, 51)
    t = 1

    F_l = []

    for n in n_range:
        print(n)
    
        HCS = C_I_HCS(n=n,
            thresholds=[c],
            censoring_numbers=[n],
            Digits=50,
            Method="analytic")
        

        F_l.append(HCS.cdf_MLE(t=t, theta=theta))


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.5))

    ax.plot(n_range, F_l, color="red", marker="o", ms=3, label=r"$F^{MLE}(t)$")

    ax.set_xlabel("n")
    ax.set_ylabel("CDFs of MLE")
    ax.legend()

    plt.tight_layout()
    plt.show()




"""Idea: Kommt die Oszillation von F^TTT()?"""
def cdf_of_TTT(t = 1.0, c = 1, theta = 1):
    p = 1 - np.exp(-c * theta)
    n_range = np.arange(4, 20 + 1)
    cdf_val_TTT = []

    for n in n_range:
        print(n)
        weights = binom.pmf(np.arange(n + 1), n, p)  # P(R = r), r = 0..n
        
        total = weights[0]  # r = 0 term: exp(-n) * 1
        HCS = C_I_HCS(
                n=n,
                thresholds=[c],
                censoring_numbers=[n],
                Digits=100,
                Method="analytic",
            )
        for r in range(1, n + 1):
            total += weights[r] * (1 - HCS.cdf_TTT(t=r/t, theta = theta))

        cdf_val_TTT.append(total)

    _, ax = plt.subplots()
    ax.plot(n_range, cdf_val_TTT, "o-")
    ax.set_xlabel("n")
    ax.set_ylabel("CDF of TTT (mixture over r)")
    ax.set_title(f"Mixture CDF of TTT vs. sample size (t = {t}, c = {c}, theta = {theta})")
    plt.show()

    return n_range, np.array(cdf_val_TTT)

if __name__ == "__main__":
    _, res = cdf_of_TTT(t = 1, c = 0.5, theta = 10)
    print(res)
    #plot_influence_of_CDF_MLE()