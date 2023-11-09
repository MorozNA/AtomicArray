import logging
import matplotlib.pyplot as plt
import numpy as np
from src.radiative_shift import DisorderedComb
from src.radiative_shift import CubicComb
from src.radiative_shift.dyson_solvers.param import LBAR


logging.basicConfig(format='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO)

period = 852e-7 / 2 / LBAR
a = period
height = 0.5 * a
width = 0.5 * a
length_etched = 1.5 * a
length = period * 5
n0 = 20
density = n0  # / (852e-7 / 2 / np.pi) ** 3

print(int(length / period))
print(round(length / period))

num_etched = int(length / period)
num_atoms = int(length * height * width * density) + int(num_etched * (length_etched * height * (period / 2) * density))
print(num_atoms)

test = DisorderedComb(length, period, density)
num_etched = int(length / period)
V = (length * height * width + num_etched * (length_etched * height * (period / 2)))
n0_after = len(test.x) / V * (852e-7 / 2 / np.pi) ** 3
print(n0_after)

test.add_atom(0.25 * 200e-7)
# test.replace_atom_z(0.25*length)

x = test.x
y = test.y
z = test.z
print(len(x), len(y), len(z))

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_aspect('auto')

scat = ax.scatter(z, x, y, s=3)
ax.set_xlabel('zlabel')
ax.set_ylabel('xlabel')
ax.set_zlabel('ylabel')
set_axes_equal(ax)
plt.show()

print(width / length)
plt.plot(test.z[:-1] / period, test.x[:-1] / period, '.')
plt.xlabel('x / a')
plt.show()