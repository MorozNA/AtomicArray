import matplotlib.pyplot as plt
from src.radiative_shift import HexagonSphere
import numpy as np

dens = 1 * (2 * 3.14 / 780e-7)**3
radius = 10/(2 * 3.14 / 780e-7)

# Tessting execution
test = HexagonSphere(dens, radius)

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
factors = np.linspace(0.9, 1.5, 100)
density = []
natoms = []
for factor in factors:
    model = HexagonSphere(dens, radius * factor)
    density.append(len(model.x) / 4 / np.pi * 3 / (radius * factor)**3)
    natoms.append(len(model.x))

plt.plot(factors, density)
plt.show()
