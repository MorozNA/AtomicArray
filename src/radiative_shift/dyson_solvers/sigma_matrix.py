import numpy as np
from .param import DDI, HBAR, GAMMA, C, OM, N_refr, KV
from src.radiative_shift.model import GeneralModel
from src.radiative_shift.tools import d_up, d_down, d_up_v, d_down_v
from src.radiative_shift.tools import find_kd, reshape_to_blocks
from abc import ABC


class SigmaMatrix(ABC):
    sigma: np.array


class MarkovianSigmaMatrixForV(SigmaMatrix):
    """
    Class for calculating self-energy part for a system of n atoms.
    """

    def __init__(self, model: GeneralModel, kd=None):
        if kd is None:
            self.kd = (OM - find_kd(N_refr, model.properties.density) * GAMMA) / C
        else:
            self.kd = kd

        xm, x0, rr = model.calculate_distances()
        nat = model.properties.noa

        x = np.zeros((3, len(xm), len(xm)), dtype='complex128')
        x[0] = xm
        x[1] = x0
        x[2] = -np.conj(xm)

        g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # Sigma V - V

        d1 = ((DDI * 1 - 1j * self.kd * rr - (self.kd * rr) ** 2) / ((rr + np.identity(nat)) ** 3) *
              np.exp(1j * self.kd * rr)) * (np.ones(nat) - np.identity(nat))
        d2 = -1 * ((DDI * 3 - 3 * 1j * self.kd * rr - (self.kd * rr) ** 2) / ((rr + np.identity(nat)) ** 3) *
                   np.exp(1j * self.kd * rr)) * (np.ones(nat) - np.identity(nat))

        di = np.zeros([nat, nat, 3, 3], dtype=np.complex)

        # Calculate the dipole elements for each value of m
        m = [-1, 0, 1]
        u = np.array([d_up_v(self.kd, 0, mi) for mi in m])
        v = np.array([d_down_v(self.kd, 0, mi) for mi in m])
        for i in range(len(m)):
            for j in range(len(m)):
                di[:, :, i, j] = np.dot(g @ u[i], v[j]) * d1
                outer = np.outer(u[i], v[j])
                for k in range(3):
                    for l in range(3):
                        di[:, :, i, j] = di[:, :, i, j] + outer[k, l] * x[k] * x[l] * d2

        self.sigma = np.reshape(np.transpose(di, axes=(0, 2, 1, 3)), [3 * nat, 3 * nat])

    def get_resolvent_for_v(self, omega):
        needed = (omega - self.kd * C + 1j * GAMMA / 2) * np.eye(len(self.sigma)) - self.sigma / HBAR
        resolvent = np.linalg.inv(needed)
        return resolvent

    def get_sigma_outside(self, model: GeneralModel, resolvent, m0, m, F0=1, F=0, J0=1 / 2, J=3 / 2, I=3 / 2):
        # TODO: add commentaries, optimize nested loops
        xm, x0, rr = model.calculate_distances_to_atom()
        xm = xm[:-1]
        x0 = x0[:-1]
        rr = rr[:-1]
        nat = model.properties.noa

        x = np.zeros((3, len(xm)), dtype='complex128')
        x[0] = xm
        x[1] = x0
        x[2] = -np.conj(xm)
        g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        d1 = ((DDI * 1 - 1j * KV * rr - (KV * rr) ** 2) / (rr ** 3) *
              np.exp(1j * KV * rr))
        d2 = -1 * ((DDI * 3 - 3 * 1j * KV * rr - (KV * rr) ** 2) / (rr ** 3) *
                   np.exp(1j * KV * rr))

        mV = [-1, 0, 1]

        dc = np.zeros([1, nat * len(m0), len(m), len(mV)], dtype=np.complex)
        db = np.zeros([1, nat * len(m0), len(mV), len(m)], dtype=np.complex)

        uc = np.zeros((len(m0), len(m), 3), dtype=np.ndarray)
        vc = np.array([d_down_v(self.kd, 0, mi) for mi in mV])

        vb = np.zeros((len(m0), len(m), 3), dtype=np.ndarray)
        ub = np.array([d_up_v(self.kd, 0, mi) for mi in mV])

        for i in range(len(m0)):
            for j in range(len(m)):
                uc[i, j] = d_up(m0[i], m[j], F0, F, J0, J, I)
                vb[i, j] = d_down(m0[i], m[j], F0, F, J0, J, I)
                for k in range(len(mV)):
                    dc[:, i * nat:(i + 1) * nat, j, k] = np.dot(g @ uc[i, j], vc[k]) * d1
                    db[:, i * nat:(i + 1) * nat, k, j] = np.dot(g @ ub[k], vb[i, j]) * d1
                    # TODO: check whether is true
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

    def get_sigma_outside_v(self, model: GeneralModel, resolvent):
        m0 = [0]
        m = [-1, 0, 1]

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

        uc = np.array([d_up_v(KV, 0, mi) for mi in mV])
        vc = np.array([d_down_v(self.kd, 0, mi) for mi in mV])

        vb = np.array([d_down_v(KV, 0, mi) for mi in mV])
        ub = np.array([d_up_v(self.kd, 0, mi) for mi in mV])

        for i in range(len(m)):
            for j in range(len(m)):
                dc[:, :, i, j] = np.dot(g @ uc[i], vc[j]) * d1
                db[:, :, i, j] = np.dot(g @ ub[i], vb[j]) * d1
                outerc = np.outer(uc[i], vc[j])
                outerb = np.outer(ub[i], vb[j])
                for k in range(3):
                    for l in range(3):
                        dc[:, :, i, j] = dc[:, :, i, j] + outerc[k, l] * x[k] * x[l] * d2
                        db[:, :, i, j] = db[:, :, i, j] + outerb[k, l] * x[k] * x[l] * d2

        sigmaC = np.reshape(np.transpose(dc, axes=(0, 2, 1, 3)), [len(m), len(m0) * nat * len(mV)])
        sigmaB = np.reshape(np.transpose(db, axes=(1, 2, 0, 3)), [len(m0) * nat * len(mV), len(m)])

        sigma_out = -1j * GAMMA / 2 * np.identity(len(m))
        for i in range(len(m0)):
            c = sigmaC[:, i * nat * len(mV): (i + 1) * nat * len(mV)] / HBAR
            b = sigmaB[i * nat * len(mV): (i + 1) * nat * len(mV), :] / HBAR
            sigma_out = sigma_out + c @ resolvent @ b
        return sigma_out

    def t_matrix(self, k1, e1, k2, e2, omega, model: GeneralModel):
        m = [-1, 0, 1]
        # TODO: check dimensions of resolvent
        resolvent = reshape_to_blocks(self.get_resolvent_for_v(omega), 3, 3) / HBAR  # get_res gives (res * HBAR)

        # TODO: write tests
        nat = model.properties.noa
        constant = 2 * np.pi * HBAR * np.sqrt(np.linalg.norm(k1) * C * np.linalg.norm(k2) * C)  # / V

        de1 = np.array([d_up_v(self.kd, 0, mi) @ e1 for mi in m])
        # TODO: check whether it is really d_up.conj()
        de2 = np.array([d_up_v(self.kd, 0, mi).conj() @ e2 for mi in m])
        outer1 = np.einsum("i,j->ij", de1, de2)
        t_ab = np.einsum("ij,abij->ab", outer1, resolvent)

        r = np.zeros((3, nat), dtype='complex128')
        r[0] = model.x
        r[1] = model.y
        r[2] = model.z
        k1r1 = np.einsum("i,ia->a", k1, r)
        k2r2 = np.einsum("i,ib->b", k2, r)
        exp1 = np.exp(1j * k1r1)
        exp2 = np.exp(-1j * k2r2)
        outer2 = np.einsum("a,b->ab", exp1, exp2)

        t = np.einsum("ab,ab", outer2, t_ab)
        return t * constant

    def total_cs(self, k1, e1, omega, model: GeneralModel):
        # This is a function of omega
        # see 'coherent control of light transport ...'
        return -2 / HBAR / C * np.imag(self.t_matrix(k1, e1, k1, e1, omega, model))

    def total_cs_sum(self, k1, e1, omega, model: GeneralModel):
        # This is a function of omega
        # see 'coherent control of light transport ...'
        e2 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        cs = np.array([np.imag(self.t_matrix(k1, e1, k1, x, omega, model)) for x in e2])
        return -2 / HBAR / C * np.sum(cs)
