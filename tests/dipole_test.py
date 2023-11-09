import numpy as np
from src.radiative_shift import d_up, d_down, d_up_v, d_down_v
from src.radiative_shift.dyson_solvers.param import KV, HBAR, GAMMA

# Define Rb 87 D-2 line transition parameters
F0 = 1
F = 0
J0 = 1 / 2
J = 3 / 2
I = 3 / 2

# Calculate dipole matrix elements for kv = KV, M0 = 0, M = -1, 0, 1 using vector spherical harmonics
dm1_v = d_up_v(KV, 0, -1)
d0_v = d_up_v(KV, 0, 0)
d1_v = d_up_v(KV, 0, 1)

# Calculate dipole matrix elements for M0 = -1, 0, 1 using Wigner-Eckart theorem
dm1 = d_up(-1, 0, F0, F, J0, J, I)
d0 = d_up(0, 0, F0, F, J0, J, I)
d1 = d_up(1, 0, F0, F, J0, J, I)

# Print dipole matrix elements
print(f"Dipole matrix element for kv = KV, M0 = 0, M = -1: {dm1_v}")
print(f"Dipole matrix element for kv = KV, M0 = 0, M = 0: {d0_v}")
print(f"Dipole matrix element for kv = KV, M0 = 0, M = 1: {d1_v}")
print(f"Dipole matrix element for M0 = -1, M = 0: {dm1}")
print(f"Dipole matrix element for M0 = 0, M = 0: {d0}")
print(f"Dipole matrix element for M0 = 1, M = 0: {d1}")

print('_________________')

# Constants used for comparison
e = 4.8032e-10  # Elementary charge in CGS units
a0 = 0.5292e-8  # Bohr radius in cm

# Calculate reduced matrix element for Rb 87 transition from Steck
rb_reduced = 4.227

# TODO: add formula
# Calculate reduced matrix element for Rb 87 transition using formula
rb_reduced_calculated = np.sqrt(3 * HBAR * GAMMA / (KV ** 3)) * np.sqrt(2 * J0 + 1) / np.sqrt(2 * J + 1) / e / a0

print(rb_reduced)
print(rb_reduced_calculated)

print(np.allclose(rb_reduced, rb_reduced_calculated, rtol=1e-03))

print('\n')
print('\n')
g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])
m = [-1, 0, 1]
u = np.array([d_up_v(KV, 0, mi) for mi in m])
v = np.array([d_down_v(KV, 0, mi) for mi in m])
print(np.shape(u))
print(u[0])
print(v[0])

outer = np.outer(u[0], v[0]) / abs(u[0, 0]) / abs(v[0, 2])
outer = outer
print(np.round(outer, 2))
print(outer[1])