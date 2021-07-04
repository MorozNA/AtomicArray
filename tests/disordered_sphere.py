import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import DisorderedSphere


layers = 5
copies = 10
n = (copies + 1) * (1 + 3 * (layers + 2) * (layers + 1))

test = DisorderedSphere(n)

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

rx, rz, rr = test.getDistances()
counter = []
# В цикле повторяю всё из метода removeDuplicates, но через матрицу расстояний
for i in range(n):
    for j in range(n):
        if i != j and rr[i][j] < 0.1:
            counter.append(i)
            continue


test.removeDuplicates()
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
