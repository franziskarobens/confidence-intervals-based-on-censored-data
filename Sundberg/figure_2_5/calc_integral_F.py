"""
Fourier inversion (Levy/Gil-Pelaez-type) formula for the CDF of a sum
of N iid truncated-exponential random variables, as used in the
Type I censored exponential coverage-probability calculations
(Sundberg-type method).

    F(z|N) = lim_{T->inf} (1/pi) * Int_{-T}^{T}
                 [sin(z t / 2) / t] * exp(-i t z / 2) * phi(t)^N  dt

    phi(t) = 1/(1 - i*theta*t) *
             (1 - exp(-C*(1 - i*theta*t)/theta)) /
             (1 - exp(-C/theta))
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt



def phi(t, theta, C):
    """Characteristic function of Exp(theta) truncated to [0, C]."""
    t = np.asarray(t, dtype=complex)
    numerator = 1 - np.exp(-C * (1 - 1j * theta * t) / theta)
    denominator = 1 - np.exp(-C / theta)
    return (1 / (1 - 1j * theta * t)) * numerator / denominator


def fourier_inversion_integral(T, z, N, theta, C):
    """
    Computes F(z|N) via the Fourier inversion formula above,
    truncating the (formally infinite) integral at +-T.

    Parameters
    ----------
    T : float
        Truncation point of the integral (should be large; the
        true formula takes T -> infinity).
    z : float
        Evaluation point of the CDF (e.g. N*t_thresh - (n-N)*C).
    N : int
        Power to which the characteristic function is raised
        (number of summed truncated-exponential variables).
    theta : float
        Mean of the (untruncated) exponential distribution.
    C : float
        Truncation point of the underlying exponential.

    Returns
    -------
    float
        The value of F(z|N), clipped to [0, 1] to guard against
        small numerical overshoot from the oscillatory integrand.
    """

    if N == 0:
        return 0.5 * np.sign(z)

    def integrand_real(t):
        # Handle the removable singularity at t = 0
        # (limit of sin(z t/2)/t as t->0 is z/2, and phi(0)^N = 1)
        if abs(t) < 1e-12:
            return z / 2.0
        val = (np.sin(z * t / 2.0) / t) * np.exp(-1j * t * z / 2.0) * phi(t, theta, C) ** N
        return val.real

    # phi(-t) = conj(phi(t)) and sin(zt/2)/t is even in t, so the
    # integrand satisfies g(-t) = conj(g(t)); hence the integral over
    # [-T, T] equals 2 * Int_0^T Re(g(t)) dt, avoiding a complex quad.
    integral, _ = integrate.quad(
        integrand_real, 0, T,
        limit=400,       # allow many subintervals for the oscillatory integrand
        epsabs=1e-10,
        epsrel=1e-10,
    )

    F = (2.0 / np.pi) * integral
    # return float(min(max(F, 0.0), 1.0))
    return float(F)


if __name__ == "__main__":
    # Example / sanity check

    actual_1 = fourier_inversion_integral(z = 0, N = 10, C = 1, theta = 1, T = 300)
    actual_2 = fourier_inversion_integral(z = 0, N = 1, C = 1, theta = 1, T = 300)
    actual_3 = fourier_inversion_integral(z = 1, N = 0, C = 1, theta = 1, T = 1)
    actual_4 = fourier_inversion_integral(z = 1, N = 0, C = 1, theta = 1, T = 2)
    actual_5 = fourier_inversion_integral(z = 1, N = 0, C = 1, theta = 1, T = 200)
    actual_6 = fourier_inversion_integral(z = 1, N = 0, C = 1, theta = 1, T = 5000)
    actual_7 = fourier_inversion_integral(z = -1, N = 0, C = 1, theta = 1, T = 5000)
    

    expected_1 = 0
    expected_2 = 0
    expected_3 = 0.946083 / np.pi
    expected_4 = 1.60541 / np.pi
    expected_5 = 1.56838 / np.pi
    expected_6 = 0.5
    expected_7 = - 0.5


    print(f"actual_1 = {actual_1:.6f}, expected_1 = {expected_1:.6f}")
    print(f"actual_2 = {actual_2:.6f}, expected_2 = {expected_2:.6f}")                
    print(f"actual_3 = {actual_3:.6f}, expected_3 = {expected_3:.6f}")
    print(f"actual_4 = {actual_4:.6f}, expected_4 = {expected_4:.6f}")
    print(f"actual_5 = {actual_5:.6f}, expected_5 = {expected_5:.6f}")
    print(f"actual_6 = {actual_6:.6f}, expected_6 = {expected_6:.6f}")
    print(f"actual_7 = {actual_7:.6f}, expected_7 = {expected_7:.6f}")


    T = 100
    N = 2
    theta = 10
    C = 20
   
    z_values = np.linspace(-1, 40, 400)
    F_values = [fourier_inversion_integral(T=T, z=z, N=N, theta=theta, C=C) for z in z_values]

    fig, ax = plt.subplots()
    ax.plot(z_values, F_values, linewidth=1.5)
    ax.set_xlabel("z")
    ax.set_ylabel("F(z|N)")
    ax.set_title(f"fourier_inversion_integral(N={N}, C={C}, theta={theta}, T={T})")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("fourier_integral_plot.png", dpi=150)