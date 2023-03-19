import matplotlib.pyplot as plt
import numpy as np
from src.radiative_shift import DisorderedModel


l = 1.0 * 2 * np.pi * (780e-7 / 2 / np.pi)
r = 200 / 780 * 2 * np.pi * (780e-7 / 2 / np.pi)
n0 = 10
density = n0 / (780e-7 / 2 / np.pi) ** 3

test = DisorderedModel(l, r, density)

x = test.x
y = test.y
z = test.z

xm, x0, xr = test.calculate_distances()
print(xm[-1])
xr[-1, -1] = 0
xr = xr ** 6
print('xr = ', xr[-1])

den_disordered = len(x) / (np.pi * l * r ** 2)
print(den_disordered / density)

print(len(x))
assert len(x) == len(y)
assert len(x) == len(z)

plt.scatter(x, y, s=3)
plt.axis('square')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=3)
plt.show()