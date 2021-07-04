import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import HexagonModel

# Тест программы.

layers = 5
copies = 10

# Перенес в конструктор, пометил addLayer и addCopy чертами

test = HexagonModel(layers, copies)

x = test.x
y = test.y
z = test.z

assert len(x) == len(y)
assert len(x) == len(z)
assert len(x) == (copies + 1) * (1 + 3 * (layers + 2) * (layers + 1))

plt.scatter(x, y, s=5)
plt.axis('square')
plt.axis('off')
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z, s=5)
plt.show()

