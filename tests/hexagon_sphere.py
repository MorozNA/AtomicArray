import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import HexagonSphere

layers = 7
copies = 11


test = HexagonSphere(layers, copies)

x = test.x
y = test.y
z = test.z

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x, y, z)
plt.show()
