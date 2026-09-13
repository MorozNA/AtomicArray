import logging
import matplotlib.pyplot as plt
import numpy as np
from src.radiative_shift import DisorderedSphere
from src.radiative_shift.atomspecies import AtomSpecies

logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=38.11e6)
LBAR = atom.lbar
density = 20 / LBAR ** 3
radius = 200e-7
n = int(density * 4 / 3 * np.pi * radius ** 3)

test = DisorderedSphere(density, radius, atom, atom)

x = test.x
y = test.y
z = test.z

print(n)
assert len(x) == len(y)
assert len(x) == len(z)
assert len(x) == n

plt.scatter(x, y, s=3)
plt.axis('square')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=3)
plt.show()

rx, rz, rr = test.calculate_distances()
counter = []
lfactor = density ** (-1/3) / 2
# Keep the first atom in each nearby pair, as remove_near_duplicates does.
for i in range(n):
    for j in range(i):
        if j not in counter and rr[i][j] < lfactor:
            counter.append(i)
            break


test.remove_near_duplicates(lfactor / LBAR)
x = test.x
y = test.y
z = test.z

# Check whether results agree
assert len(x) == n - len(set(counter))
_, _, rr = test.calculate_distances()
assert np.all(rr[np.triu_indices(len(x), k=1)] >= lfactor)

plt.scatter(x, y, s=3)
plt.axis('square')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=3)
plt.show()
