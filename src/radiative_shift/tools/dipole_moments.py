import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
# from sympy.physics.quantum.cg import CG
from src.radiative_shift.dyson_solvers.param import HBAR, GAMMA, KV


# 133Cs parameters are F0=4, F=5, J0=1/2, J=3/2, I=7/2
# 87Rb parameters are F0=1, F=0, J0=1/2, J=3/2, I=3/2
def d_up(M0, M, F0=4, F=5, J0=1 / 2, J=3 / 2, I=7 / 2):
    """
    Calculates the dipole matrix element <J0,M0|q|J,M> using the Wigner-Eckart theorem.

    Arguments:
    M0 (float): The magnetic quantum number for the ground state.
    M (float): The magnetic quantum number for the excited state.
    F0 (float, optional): The hyperfine quantum number for the ground state.
    F (float, optional): The hyperfine quantum number for the excited state.
    J0 (float, optional): The total angular momentum for the ground state.
    J (float, optional): The total angular momentum for the excited state.
    I (float, optional): The nuclear spin.

    Returns:
    numpy.ndarray: The contravariant spherical components of the dipole operator for chosen transition.
    """
    q = [-1, 0, 1]
    j3 = np.array([wigner_3j(F, 1, F0, M, q[i], -M0) for i in range(3)])
    j6 = wigner_6j(J0, J, 1, F, F0, I)
    d_vec = (-1) ** (J0 + M0 + I) * np.sqrt((2 * F0 + 1) * (2 * F + 1)) * j6 * j3
    reduced_up = np.sqrt(3 * HBAR * GAMMA / 4 / (KV ** 3)) * np.sqrt(2 * J + 1)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex) * reduced_up


def d_down(M0, M, F0=4, F=5, J0=1 / 2, J=3 / 2, I=7 / 2):
    d_vec = d_up(M0, M, F0, F, J0, J, I)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex)


def d_up_v(k, M0, M):
    J0, J = 0, 1
    q = [-1, 0, 1]
    j3 = np.array([wigner_3j(J, 1, J0, M, q[i], -M0) for i in range(3)], dtype=np.complex)
    d_vec = (-1) ** (J - 1 + M0) * j3
    reduced_v = np.sqrt(3 * HBAR * GAMMA / 4 / (k ** 3)) * np.sqrt(2 * J + 1)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex) * reduced_v


def d_down_v(k, M0, M):
    d_vec = d_up_v(k, M0, M)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex)
