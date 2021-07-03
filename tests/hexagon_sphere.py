import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import HexagonSphere
from src.radiative_shift.model.param import SPHERE_RADIUS
import numpy as np

layers = 7
copies = 11

# Tessting execution
test = HexagonSphere(layers, copies)

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
RADIUS = SPHERE_RADIUS
factors = np.linspace(0.9, 1.5, 100)
density = []
natoms = []
for factor in factors:
    model = HexagonSphere(layers, copies, RADIUS * factor)
    density.append(len(model.x) / 4 / np.pi * 3 / (RADIUS * factor)**3)
    natoms.append(len(model.x))

plt.plot(factors, density)
plt.show()
