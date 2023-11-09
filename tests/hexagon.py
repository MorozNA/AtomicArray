import logging
import numpy as np
import matplotlib.pyplot as plt
from src.radiative_shift import HexagonModel
from src.radiative_shift.dyson_solvers.param import LBAR

# Тест программы.

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

l = 2 * 780e-7 / LBAR
r = 200e-7 / LBAR
density = 15  # * LBAR ** 3
n = int(density * l * np.pi * r ** 2)
print(n)

test = HexagonModel(l, r, density)
test.rotate_model(np.pi/6)

test.add_atom(1.2 * r)

x = test.x
y = test.y
z = test.z

assert len(x) == len(y)
assert len(x) == len(z)

plt.scatter(x, y, s=5)
plt.axis('square')
plt.axis('off')

circle1 = plt.Circle((0, 0), r, color='r', fill=False)
plt.gca().add_patch(circle1)

# plt.xlim(-1.3 * test.properties.radius, 1.3 * test.properties.radius)
# plt.ylim(-1.3 * test.properties.radius, 1.3 * test.properties.radius)
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=5)
ax.scatter3D(x[-1], y[-1], z[-1], s=15)
plt.show()
