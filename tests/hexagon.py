import logging

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import HexagonModel

# Тест программы.

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

n = 2400
r = 4.8

test = HexagonModel(n, r)

x = test.x
y = test.y
z = test.z

assert len(x) == len(y)
assert len(x) == len(z)

plt.scatter(x, y, s=5)
plt.axis('square')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=5)
plt.show()
