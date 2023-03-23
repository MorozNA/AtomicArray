import numpy as np
from .param import KV, DDI, HBAR, GAMMA, C
from src.radiative_shift.model import GeneralModel
from src.radiative_shift.dyson_solvers.dipole_moments import d_up, d_down, d_up_v, d_down_v
from abc import ABC


class SigmaMatrix(ABC):
    sigma: np.array


class MarkovianSigmaMatrixForV(SigmaMatrix):
    """
    Class for calculating self-energy part for a system of n atoms.
    """
    def __init__(self, model: GeneralModel, kd):
        xm, x0, rr = model.calculate_distances()
        nat = np.shape(xm)[0]

        x = [xm, x0, -np.conj(xm)]
        g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # Sigma V - V

        d1 = ((DDI * 1 - 1j * KV * rr - (KV * rr) ** 2) / ((rr + np.identity(nat)) ** 3) *
              np.exp(1j * KV * rr)) * (np.ones(nat) - np.identity(nat))
        d2 = -1 * ((DDI * 3 - 3 * 1j * KV * rr - (KV * rr) ** 2) / ((rr + np.identity(nat)) ** 3) *
                   np.exp(1j * KV * rr)) * (np.ones(nat) - np.identity(nat))

        di = np.zeros([nat, nat, 3, 3], dtype=np.complex)

        # Calculate the dipole elements for each value of m
        m = [-1, 0, 1]
        u = np.array([d_up_v(kd, 0, mi) for mi in m])
        v = np.array([d_down_v(kd, 0, mi) for mi in m])
        for i in range(len(m)):
            for j in range(len(m)):
                di[:, :, i, j] = np.dot(g @ u[i], v[j]) * d1
                outer = np.outer(u[i], v[j])
                for k in range(3):
                    for l in range(3):
                        di[:, :, i, j] = di[:, :, i, j] + outer[k, l] * x[k] * x[l] * d2

        self.sigmaV = np.reshape(np.transpose(di, axes=(0, 2, 1, 3)), [3 * nat, 3 * nat])

    def get_resolvent_for_v(self, omega, kd):
        needed = (omega - kd * C + 1j * GAMMA / 2) * np.eye(len(self.sigmaV)) - self.sigmaV / HBAR
        resolvent = np.linalg.inv(needed)
        return resolvent

    def get_simga_outside(self, resolvent, m0, m):
        pass