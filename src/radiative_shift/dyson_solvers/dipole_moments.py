import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from sympy.physics.quantum.cg import CG
from src.radiative_shift.dyson_solvers.param import HBAR, GAMMA, KV

# 133Cs parameters are F0=4, F=5, J0=1/2, J=3/2, I=7/2
# 87Rb parameters are F0=1, F=0, J0=1/2, J=3/2, I=3/2
def d_up(M0, M, F0=4, F=5, J0=1/2, J=3/2, I=7/2):
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
    numpy.ndarray: The Contravariant spherical components of the dipole operator for chosen transition.
    """
    q = [-1, 0, 1]
    j6 = wigner_6j(J, F, I, F0, J0, 1)
    j3 = np.array([wigner_3j(F0, 1, F, M0, q[i], -M) for i in range(3)])
    d_vec = (-1) ** (I + J - M) * np.sqrt((2 * F0 + 1) * (2 * F + 1)) * j6 * j3
    d = np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex)
    reduced_up = np.sqrt(3 * HBAR * GAMMA / (KV ** 3))
    return d * reduced_up


def d_down(M0, M, F0=4, F=5, J0=1 / 2, J=3 / 2, I=7 / 2):
    d_vec = d_up(M0, M, F0, F, J0, J, I)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex)


def d_up_v(k, M0, M):
    F0, F = 0, 1
    q = [-1, 0, 1]
    d_vec = [CG(F0, M0, 1, q[i], F, M).doit() / np.sqrt((2 * F + 1)) for i in range(3)]
    reduced_v = np.sqrt(3 * HBAR * GAMMA / (k ** 3)) * np.sqrt((2 * F + 1) / 4)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex) * reduced_v


def d_down_v(k, M0, M):
    d_vec = d_up_v(k, M0, M)
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]], dtype=np.complex)
