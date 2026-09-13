import logging
import numpy as np
import matplotlib.pyplot as plt
from src.radiative_shift import HexagonModel
from src.radiative_shift.atomspecies import AtomSpecies


logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=38.11e6)
l = 2 * 780e-7
r = 200e-7
density = 15 / atom.lbar ** 3
n = int(density * l * np.pi * r ** 2)
print(n)

test = HexagonModel(l, r, density, atom, atom)
test.rotate_about_z(np.pi/6)

test.add_medium_atom_cylindrical(1.2 * r)

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

plt.xlim(-1.3 * r, 1.3 * r)
plt.ylim(-1.3 * r, 1.3 * r)
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=5)
ax.scatter3D(x[-1], y[-1], z[-1], s=15)
plt.show()
