import numpy as np
from src.radiative_shift import EmptyModel
from src.radiative_shift import MarkovianSigmaMatrixForV
from src.radiative_shift.dyson_solvers.param import LBAR, C, GAMMA, KV


# TODO: Check whether anything changes if Gamma(V) = Gamma(\Lambda) / 3

def calc_tsc(x1, x2):
    model = EmptyModel()

    model.add_atom_xyz(x1[0], x1[1], x1[2])
    model.add_atom_xyz(x2[0], x2[1], x2[2])
    # TODO: resolve problem with measure_properties() after atom addition of new atoms
    model.measure_properties()

    sigma_v = MarkovianSigmaMatrixForV(model, KV)

    x = np.linspace(-15, 15, 1000)
    tcs = []
    e1 = np.array([1, 0, 0])
    for i in x:
        omega = sigma_v.kd * C + i * GAMMA
        k1 = omega / C * np.array([0, 0, 1])
        tcs.append(sigma_v.total_cs(k1, e1, omega, model) / (LBAR ** 2))
    return tcs
