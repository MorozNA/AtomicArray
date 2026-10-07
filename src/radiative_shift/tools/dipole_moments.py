import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import HBAR


# 133Cs parameters are F0=4, F=5, J0=1/2, J=3/2, I=7/2
# 87Rb parameters are F0=1, F=0, J0=1/2, J=3/2, I=3/2


def dipole_mn(atom: AtomSpecies, M0, M, k=None):
    """
    Return contravariant spherical components of d_mn = <m|d|n>.

    Component order is (-1, 0, +1);
    dipole units are Gaussian CGS;
    m = |F0,M0> and n = |F,M> are the ground and excited states;
    F0 and F are ground/excited hyperfine angular momenta;
    J0 and J are ground/excited electronic angular momenta;
    I is the nuclear spin;
    M0 and M are magnetic quantum numbers.
    k is atomic wavenumber
    """

    if k is None:
        k = atom.wavenumber
    F0, F, J0, J, I = atom.F0, atom.F, atom.J0, atom.J, atom.I
    j3 = np.array([wigner_3j(F, 1, F0, M, q, -M0) for q in (-1, 0, 1)], dtype=complex)
    j6 = float(wigner_6j(J0, J, 1, F, F0, I))

    phase = (-1.0) ** (2 * F + J0 + M0 + I)
    d_vec = phase * np.sqrt((2 * F0 + 1) * (2 * F + 1)) * j6 * j3

    reduced_dipole = np.sqrt(3 * HBAR * atom.gamma * (2 * J + 1) / (4 * k**3))

    # Raise the spherical index: d^q = (-1)^q d_{-q}.
    return np.array([-d_vec[2], d_vec[1], -d_vec[0]]) * reduced_dipole


def dipole_nm(atom: AtomSpecies, M0, M, k=None):
    """Return d_nm = <n|d|m> in the same contravariant spherical convention."""
    dipole_components = dipole_mn(atom, M0, M, k)
    return np.array([-dipole_components[2], dipole_components[1], -dipole_components[0]]).conj()
