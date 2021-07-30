import numpy as np
from .param import KV, DDI, HBAR, D01M, D1M0, D00, D10, D01
from src.radiative_shift.model import GeneralModel

class markovianSigmaMatrixForV:
    """
    TODO: Напишу юниттест

    :param model:
    :return: sigma matrix
    """
    def __init__(self, model: GeneralModel):
        xm, x0, rr = model.getDistances()
        nat = np.shape(xm)[0]
        # форматирование ужасное, знаю :)
        d1 = 1 * KV * ((DDI * 1 - 1j * rr * KV - (rr * KV) ** 2) /
                       ((rr * KV + np.identity(nat)) ** 3) *
                       np.exp(1j * rr * KV)) * (np.ones(nat) - np.identity(nat))
        d2 = -1 * KV * ((DDI * 3 - 3 * 1j * rr * KV - (rr * KV) ** 2) /
                        ((KV * rr + np.identity(nat)) ** 3) * np.exp(1j * KV * rr)) * (np.ones(nat) - np.identity(nat))

        # Амплитуда перехода из |F = 0; F = 1, m> в |F = 1, m'; F = 0>
        di = np.zeros([nat, nat, 3, 3], dtype=np.complex)

        di[:, :, 0, 0] = D01M * D1M0 * (d1 + np.conj(xm) * xm * d2)
        di[:, :, 0, 1] = D01M * D00 * (np.conj(xm) * x0 * d2)
        di[:, :, 0, 2] = D01M * D10 * (-np.conj(xm) * np.conj(xm) * d2)
        di[:, :, 1, 0] = D00 * D1M0 * (x0 * xm * d2)
        di[:, :, 1, 1] = D00 * D00 * (d1 + x0 * x0 * d2)
        di[:, :, 1, 2] = D00 * D10 * (-x0 * np.conj(xm) * d2)
        di[:, :, 2, 0] = D01 * D1M0 * (-xm * xm * d2)
        di[:, :, 2, 1] = D01 * D00 * (-xm * x0 * d2)
        di[:, :, 2, 2] = D01 * D10 * (d1 + xm * np.conj(xm) * d2)

        self.sigma = np.zeros([3 * nat, 3 * nat], dtype=np.complex)
        #  Тоже лучше через numpy
        for n1 in range(nat):  # initial excited
            for n2 in range(nat):  # final excited
                for i in range(3):
                    for j in range(3):
                        self.sigma[3 * n1 + j, 3 * n2 + i] = 3 * di[n1, n2, j, i]

    def heWhoShallNotBeNamed(self, chi):
        needed = np.linalg.inv((np.add(np.eye(len(self.sigma)), chi * self.sigma)))
        return needed
