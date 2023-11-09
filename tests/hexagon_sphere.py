import logging
import matplotlib.pyplot as plt
from src.radiative_shift import HexagonSphere
import numpy as np
from src.radiative_shift.dyson_solvers.param import LBAR, KV

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

R = 250
DEN = 20

r = R / 780 * 2 * np.pi * LBAR
density = DEN * KV ** 3
test = HexagonSphere(r, density)

x = test.x
y = test.y
z = test.z

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z)
plt.show()

# ax, ay, az = test.getDistances()  # Возвращает 4 матрицы np.array
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
