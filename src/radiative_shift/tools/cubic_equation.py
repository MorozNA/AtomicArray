import numpy as np
from src.radiative_shift.dyson_solvers.param import OM, GAMMA, N_refr


def solve_qubic(a, b, c, d):
    # https://math.stackexchange.com/questions/18867/how-does-one-solve-a-cubic-polynomial-with-complex-coefficients
    # https://math.stackexchange.com/questions/15865/why-not-write-the-solutions-of-a-cubic-this-way/18873#18873
    p = b / a
    q = c / a
    f = d / a

    A1 = np.sqrt(- p ** 2 * q ** 2 + 4 * q ** 3 + 4 * p ** 3 * f - 18 * p * q * f + 27 * f ** 2)
    A = (- 2 * p ** 3 + 9 * p * q - 27 * f + 3 * np.sqrt(3) * A1) ** (1 / 3) / (3 * 2 ** (1 / 3))
    B = (-p ** 2 + 3 * q) / 9 / A

    x1 = -p / 3 + A - B
    x2 = - p / 3 + (-1 - 1j * np.sqrt(3)) / 2 * A - (-1 + 1j * np.sqrt(3)) / 2 * B
    x3 = - p / 3 + (-1 + 1j * np.sqrt(3)) / 2 * A - (-1 - 1j * np.sqrt(3)) / 2 * B
    return x1, x2, x3


def find_kd(n_refr, n0):
    # ! here n0 == n0 * LBAR ** 3
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

        x = solve_qubic(a, b, c, d)
        eps.append(x[0] ** 2)

    # TODO: добавить assert близко ли к нулю попали
    to_solve = abs(np.real(np.array(eps)) - (n_refr ** 2))
    k_del = del_omega[np.argmin(to_solve)]
    return k_del


def find_eps(n_refr, n0, omega):
    # ! here n0 == n0 * LBAR ** 3
    F = 1
    F0 = 0

    alpha = 1  # (2 * F0 + 1) / 3
    ro = n0 * (2 * F + 1) / 3 / (2 * F0 + 1)

    # TODO: check this
    om_d = (OM - find_kd(N_refr, n0) * GAMMA)
    del_omega = (omega - om_d) / GAMMA

    a = 1j / 2
    b = alpha * np.pi * ro + del_omega
    c = - 1j / 2
    d = alpha * 2 * np.pi * ro - del_omega

    x = solve_qubic(a, b, c, d)
    eps = x[0] ** 2
    return eps
