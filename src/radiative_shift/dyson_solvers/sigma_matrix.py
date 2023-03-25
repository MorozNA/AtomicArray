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

    def get_sigma_outside(self, model, resolvent, kd, m0, m, F0=1, F=0, J0=1/2, J=3/2, I=3/2):
        # TODO: add commentaries, optimize nested loops
        xm, x0, rr = model.calculate_distances_to_atom()
        xm = xm[:-1]
        x0 = x0[:-1]
        rr = rr[:-1]

        x = [xm, x0, -np.conj(xm)]
        g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        nat = len(xm)

        d1 = ((DDI * 1 - 1j * KV * rr - (KV * rr) ** 2) / (rr ** 3) *
              np.exp(1j * KV * rr))
        d2 = -1 * ((DDI * 3 - 3 * 1j * KV * rr - (KV * rr) ** 2) / (rr ** 3) *
                   np.exp(1j * KV * rr))

        mV = [-1, 0, 1]

        dc = np.zeros([1, nat * len(m0), len(m), len(mV)], dtype=np.complex)
        db = np.zeros([1, nat * len(m0), len(mV), len(m)], dtype=np.complex)

        uc = np.zeros((len(m0), len(m), 3), dtype=np.ndarray)
        vc = np.array([d_down_v(kd, 0, mi) for mi in mV])

        vb = np.zeros((len(m0), len(m), 3), dtype=np.ndarray)
        ub = np.array([d_up_v(kd, 0, mi) for mi in mV])

        for i in range(len(m0)):
            for j in range(len(m)):
                uc[i, j] = d_up(m0[i], m[j], F0, F, J0, J, I)
                vb[i, j] = d_down(m0[i], m[j], F0, F, J0, J, I)
                for k in range(len(mV)):
                    dc[:, i * nat:(i + 1) * nat, j, k] = np.dot(g @ uc[i, j], vc[k]) * d1
                    db[:, i * nat:(i + 1) * nat, k, j] = np.dot(g @ ub[k], vb[i, j]) * d1
                    outerc = np.outer(uc[i, j], vc[k])
                    outerb = np.outer(ub[k], vb[i, j])
                    for q in range(3):
                        for p in range(3):
                            dc[:, i * nat:(i + 1) * nat, j, k] = \
                                dc[:, i * nat:(i + 1) * nat, j, k] + outerc[q, p] * x[q] * x[p] * d2
                            db[:, i * nat:(i + 1) * nat, k, j] = \
                                db[:, i * nat:(i + 1) * nat, k, j] + outerb[q, p] * x[q] * x[p] * d2

        sigmaC = np.reshape(np.transpose(dc, axes=(0, 2, 1, 3)), [len(m), len(m0) * nat * len(mV)])
        sigmaB = np.reshape(np.transpose(db, axes=(1, 2, 0, 3)), [len(m0) * nat * len(mV), len(m)])

        sigma_out = -1j * GAMMA / 2 * np.identity(len(m))
        for i in range(len(m0)):
            c = sigmaC[:, i * nat * len(mV): (i + 1) * nat * len(mV)] / HBAR
            b = sigmaB[i * nat * len(mV): (i + 1) * nat * len(mV), :] / HBAR
            sigma_out = sigma_out + c @ resolvent @ b
        return sigma_out
