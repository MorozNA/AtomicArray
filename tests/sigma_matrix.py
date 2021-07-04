from src.radiative_shift import markovianSigmaMatrixForV
from src.radiative_shift import HexagonSphere
from src.radiative_shift.model.param import SPHERE_RADIUS as RADIUS
import numpy as np


x = np.linspace(0.1, 1.5, 10)
y = []
for i in np.linspace(0.1, 1.5, 20):
    radius = RADIUS * i
    model = HexagonSphere(10, 10, radius)
    natoms = len(model.x)
    dens = 3 * natoms / (4 * np.pi * radius**3)
    print(dens)

    eigs, _ = np.linalg.eig(markovianSigmaMatrixForV(model))
    eigs = np.sort(np.imag(eigs))
    print("0 = {:1.5}".format(abs(eigs[-1] - 0.5))) # dark state exists
    print("1 = {:2.2}".format(-1*(radius**2 * eigs[0] / natoms))) # see PRA 93, 043830
    y.append(-1*(radius**2 * eigs[0] / natoms))

from matplotlib import pyplot as plt
plt.plot(x,y)
plt.show()
