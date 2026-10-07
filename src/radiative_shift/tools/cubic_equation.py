import numpy as np
from src.radiative_shift.atomspecies import AtomSpecies


def solve_cubic(a, b, c, d):
    """Return the roots of the cubic equation ax^3 + bx^2 + cx + d = 0."""
    # https://math.stackexchange.com/questions/18867/how-does-one-solve-a-cubic-polynomial-with-complex-coefficients
    # https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
    p = b / a
    q = c / a
    f = d / a

    A1 = np.sqrt(complex(- p ** 2 * q ** 2 + 4 * q ** 3 + 4 * p ** 3 * f - 18 * p * q * f + 27 * f ** 2))
    A = (- 2 * p ** 3 + 9 * p * q - 27 * f + 3 * np.sqrt(3) * A1) ** (1 / 3) / (3 * 2 ** (1 / 3))
    B = (-p ** 2 + 3 * q) / 9 / A

    x1 = -p / 3 + A - B
    x2 = - p / 3 + (-1 - 1j * np.sqrt(3)) / 2 * A - (-1 + 1j * np.sqrt(3)) / 2 * B
    x3 = - p / 3 + (-1 + 1j * np.sqrt(3)) / 2 * A - (-1 - 1j * np.sqrt(3)) / 2 * B
    return x1, x2, x3


def find_reference_detuning(n_refr, n0, atol=0.01):
    """Find (omega_0 - omega_M) / Gamma_e for the Appendix A permittivity.

    This detuning equals -delta_M / Gamma_e
    n_refr is the target refractive index; n0 is the dimensionless density
    n0 * lambda0_bar**3, scaled by the reference reduced wavelength cubed
    """
    if not (np.isfinite(n_refr) and n_refr > 1 and np.isfinite(n0) and n0 > 0):
        raise ValueError("This detuning search requires n_refr > 1 and n0 > 0")
    F = 1
    F0 = 0

    alpha = 1  # (2 * F0 + 1) / 3
    ro = n0 * (2 * F + 1) / 3 / (2 * F0 + 1)

    pr_d = np.pi * ro * (2 + n_refr ** 2) / (1 - n_refr ** 2)

    del_omega = np.linspace(pr_d - ro, pr_d + ro, 5000)  # del to gamma
    eps = []

    for k in range(len(del_omega)):
        delt = del_omega[k]

        a = 1j / 2
        b = alpha * np.pi * ro + delt
        c = - 1j / 2
        d = alpha * 2 * np.pi * ro - delt

        x = solve_cubic(a, b, c, d)
        eps.append(x[0] ** 2)

    to_solve = abs(np.real(np.array(eps)) - (n_refr ** 2))
    index = np.argmin(to_solve)

    if not np.isfinite(to_solve[index]) or to_solve[index] > atol:
        raise RuntimeError(
            f"Target permittivity was not reached: "
            f"error = {to_solve[index]:.6g}, tolerance = {atol:.6g}"
        )

    return del_omega[index]


def medium_permittivity(n_refr, n0, omega, atom_reference: AtomSpecies):
    """Return epsilon(omega) for the artificial medium of Appendix A.

    n_refr is the target refractive index; n0 is the dimensionless density
    n0 * lambda0_bar**3, scaled by the reference reduced wavelength cubed.
    The frequency offset and detuning use atom_reference.gamma, corresponding
    to the paper's choice Gamma_e = Gamma_infinity.
    """
    # ! here n0 == n0 * LBAR ** 3
    F = 1
    F0 = 0

    alpha = 1  # (2 * F0 + 1) / 3
    ro = n0 * (2 * F + 1) / 3 / (2 * F0 + 1)

    omega_M = (atom_reference.omega - find_reference_detuning(n_refr, n0) * atom_reference.gamma)
    del_omega = (omega - omega_M) / atom_reference.gamma

    a = 1j / 2
    b = alpha * np.pi * ro + del_omega
    c = - 1j / 2
    d = alpha * 2 * np.pi * ro - del_omega

    sqrt_epsilon_roots = solve_cubic(a, b, c, d)
    epsilon = sqrt_epsilon_roots[0] ** 2
    return epsilon
