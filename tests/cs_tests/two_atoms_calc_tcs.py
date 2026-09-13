import numpy as np
from src.radiative_shift import EmptyModel
from src.radiative_shift import MarkovianSigmaMatrixForV
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import C
from src.radiative_shift.tools import d_up


atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=38.11e6)
LBAR, GAMMA = atom.lbar, atom.gamma


def calc_tsc(x1, x2):
    model = EmptyModel(atom, atom)

    model.add_atom_xyz(x1[0], x1[1], x1[2])
    model.add_atom_xyz(x2[0], x2[1], x2[2])

    sigma_v = MarkovianSigmaMatrixForV(model)

    x = np.linspace(-15, 15, 1000)
    tcs = []
    e1 = np.array([1, 0, -1]) / np.sqrt(2)  # x polarization in spherical components
    de1 = np.array([np.vdot(d_up(atom, 0, mi), e1) for mi in atom.m])
    positions = np.column_stack((model.x, model.y, model.z))
    for i in x:
        omega = atom.omega + i * GAMMA
        k1 = omega / C * np.array([0, 0, 1])
        drive = np.outer(np.exp(1j * (positions @ k1)), de1).ravel()
        resolvent = sigma_v.get_resolvent_for_v(omega)
        # Optical theorem; the current resolvent already has inverse-energy units.
        tcs.append(-4 * np.pi * np.linalg.norm(k1)
                   * np.imag(np.vdot(drive, resolvent @ drive)) / LBAR ** 2)
    return tcs
