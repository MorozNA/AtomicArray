import logging
import matplotlib.pyplot as plt
from src.radiative_shift import HexagonSphere
import numpy as np
from src.radiative_shift.atomspecies import AtomSpecies

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

R = 250
DEN = 20

atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=38.11e6)
LBAR, KV = atom.lbar, atom.wavenumber
r = R / 780 * 2 * np.pi * LBAR
density = DEN * KV ** 3
test = HexagonSphere(r, density, atom, atom)

x = test.x
y = test.y
z = test.z

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z)
plt.show()

# ax, ay, az = test.calculate_distances()  # Returns 3 numpy arrays.
# print(ax)

# Testing density
# factors = np.linspace(0.9, 1.5, 100)
# density = []
# natoms = []
# for factor in factors:
#     model = HexagonSphere(radius * factor, density)
#     density.append(len(model.x) / 4 / np.pi * 3 / (radius * factor)**3)
#     natoms.append(len(model.x))
#
# plt.plot(factors, density)
# plt.show()
