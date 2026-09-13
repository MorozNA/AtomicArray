import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from src.radiative_shift import DisorderedComb
from src.radiative_shift import DisorderedModel
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import HBAR

reference_atom = AtomSpecies(F0=2, F=3, J0=1/2, J=3/2, I=3/2, lambda_nm=780, gamma=38.11e6)
medium_atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=reference_atom.gamma)
LBAR = reference_atom.lbar

period = 2 * np.pi * LBAR / 2
print(np.round(period * 10**7, 3))
a = period
height = 0.5 * a
width = 0.5 * a
length_etched = 1.5 * a
length = period * 5

num_etched = int(length / period)
V = (length * height * width + num_etched * (length_etched * height * (period / 2)))
n0 = 50
density = n0 / LBAR ** 3

num = 5000
density = num / V
print(int(V * density))

V = (length * height * width + num_etched * (length_etched * height * (period / 2)))

x = np.linspace(0.1, 5.0, 200)
# x = x[11::]
array = []

test = DisorderedComb(length, period, density, medium_atom, reference_atom)
r = 200e-7
test.set_reference_position_cylindrical(reference_atom, r*1.01)

###
# l = 2.0 * 2 * np.pi * (852e-7 / 2 / np.pi)
# r = 200 / 852 * 2 * np.pi * (852e-7 / 2 / np.pi)
# n0 = 20
# density = n0 / (852e-7 / 2 / np.pi) ** 3
# test = DisorderedModel(l, r, density, medium_atom, reference_atom)
# test.set_reference_position_cylindrical(reference_atom, r*1.01)
# V = np.pi * r ** 2 * l
# x = np.linspace(1.5, 5.0, 50)
# array = []
# print(len(test.x))
###

for i in range(len(x)):
    radius = r * x[i]
    test.set_reference_position_cylindrical(reference_atom, radius)

    _, _, xr = test.calculate_distances_to_reference_atom()

    xr = xr ** 6
    integral = np.sum(1 / xr)
    integral = V / len(xr) * integral  # The reference atom is not a medium sample.
    array.append(integral)

print(len(xr))
array = np.array(array)

# eps = 2.1025
eps = 30  # ~ \infty for gold
a0 = 0.529e-8
GAMMA = reference_atom.gamma

e = 4.8e-10
# var_rb = 36.18 * (a0 ** 2)
var_rb = 20.7 * 100 * (a0 ** 2)  # \sim sodium r^2 for n = 10
C = - 3 / (4 * np.pi * HBAR) * (eps - 1) / (eps + 2) * var_rb * e ** 2
array = C * array / GAMMA

GAMMA = 38.11
array = array * GAMMA  # To plot in MHz

plt.plot(x, array)
plt.show()

np.savetxt('./monte-carlo_rb.txt', array, fmt='%f')

# data1 = [x, array]
# data1 = np.array(data1).T
# np.savetxt('vdWRb.csv', data1, delimiter=',')
