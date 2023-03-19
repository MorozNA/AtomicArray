import matplotlib.pyplot as plt
import numpy as np
from src.radiative_shift import DisorderedSphere
from src.radiative_shift.dyson_solvers.param import LBAR


density = 100 / (LBAR ** 3)
radius = 1 * LBAR
n = int(density * 4 / 3 * np.pi * radius ** 3)

test = DisorderedSphere(density, radius)

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
# В цикле повторяю всё из метода removeDuplicates, но через матрицу расстояний
for i in range(n):
    for j in range(n):
        if i != j and rr[i][j] < lfactor:
            counter.append(i)
            continue


test.remove_duplicates(lfactor)
x = test.x
y = test.y
z = test.z

# Здесь проверяю совпадают ли результаты удаления точек
assert len(x) == n - len(set(counter))

plt.scatter(x, y, s=3)
plt.axis('square')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=3)
plt.show()
