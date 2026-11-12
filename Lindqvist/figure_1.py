import matplotlib.pyplot as plt
import numpy as np

def compute_tau(theta, U, c):
    numerator = np.sum(U / theta < c)
    denominator = np.sum(np.minimum(U / theta, c))
    return numerator / denominator



def plot_figure_1(U, c, theta_range=(0.01, 10), n_points=500):
    thetas = np.linspace(theta_range[0], theta_range[1], n_points)
    tau_values = np.array([compute_tau(th, U, c) for th in thetas])

    plt.figure()
    plt.plot(thetas, tau_values, color="blue", linewidth=2)
    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$\tau(\mathbf{U}, \theta)$")
    for u in U / c:
        plt.axvline(u, color="gray", linestyle=":")
    plt.show()


rng = np.random.default_rng() # set seed
theta = 2
n = 10
U = rng.exponential(scale=1.0, size=n)
c = 1



plot_figure_1(U=U, c=c, theta_range=(0.01, 1.3))

