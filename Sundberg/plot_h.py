# Create a helping plot for understanding the 
# theory behind the construction of C5

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

K = 0.5

x = np.linspace(1e-4, 5, 1000)
y = x - 1 - np.log(x)

# Schnittpunkte von h(x) = K exakt über Lambert W:
# x - 1 - log x = K  <=>  -x e^{-x} = -e^{-(1+K)}  <=>  x = -W(-e^{-(1+K)})
arg = -np.exp(-(1 + K))
x1 = -lambertw(arg, k=0).real    # linker Schnittpunkt  (x < 1)
x2 = -lambertw(arg, k=-1).real   # rechter Schnittpunkt (x > 1)

fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(x, y, color="C0", linewidth=2)
ax.axhline(K, color="C1", linestyle="--", linewidth=1.5)

# Vertikale Linien von der x-Achse bis zum Schnittpunkt
ax.vlines([x1, x2], 0, K, color="C3", linestyle="--", linewidth=1.5)
ax.plot([x1, x2], [K, K], "o", color="C3", zorder=3)

# Achsen: x-Achse direkt auf y = 0, y-Achse direkt auf x = 0
ax.set_xlim(0, 5)
ax.set_ylim(0, 2)
ax.margins(0)

ax.set_xticks([0, x1, 1, x2, 5])
ax.set_xticklabels(["0", r"$c_l$", "1", r"$c_u$", "5"])

ax.set_yticks([0, K, 2])
ax.set_yticklabels(["0", "K", "2"])

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$h(x)$")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("plot_h_for_C5.pdf", bbox_inches="tight")
plt.show()