import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import DisorderedModel

# Тест программы. Можно сначала добавлять копии, но тогда addLayer будет проходить по всем копиям, что занимает время.

n = 1400
r = 6
density = 40

# Перенес в конструктор, пометил addLayer и addCopy чертами

test = DisorderedModel(n, r, density)

x = test.x
y = test.y
z = test.z

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
