import numpy as np
from src.radiative_shift import dipole_mn, dipole_nm
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import HBAR

# Define Rb 87 D-2 line transition parameters
F0 = 1
F = 0
J0 = 1 / 2
J = 3 / 2
I = 3 / 2

atom = AtomSpecies(F0=F0, F=F, J0=J0, J=J, I=I, lambda_nm=780, gamma=38.11e6)
v_atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=atom.gamma)
KV, GAMMA = atom.wavenumber, atom.gamma

# The general dipole functions also handle the V transition.
dm1_v = dipole_mn(v_atom, 0, -1)
d0_v = dipole_mn(v_atom, 0, 0)
d1_v = dipole_mn(v_atom, 0, 1)

# Calculate dipole matrix elements for M0 = -1, 0, 1 using Wigner-Eckart theorem
dm1 = dipole_mn(atom, -1, 0)
d0 = dipole_mn(atom, 0, 0)
d1 = dipole_mn(atom, 1, 0)

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

# Calculate reduced matrix element for Rb 87 transition using formula
rb_reduced_calculated = np.sqrt(3 * HBAR * GAMMA * (2 * J + 1) / (4 * KV ** 3)) / np.sqrt(2 * J0 + 1) / e / a0

print(rb_reduced)
print(rb_reduced_calculated)

print(np.allclose(rb_reduced, rb_reduced_calculated, rtol=1e-03))
np.testing.assert_allclose(rb_reduced, rb_reduced_calculated, rtol=1e-03)

print('\n')
print('\n')
g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])
m = [-1, 0, 1]
u = np.array([dipole_mn(v_atom, 0, mi) for mi in m])
v = np.array([dipole_nm(v_atom, 0, mi) for mi in m])
d_unit = np.sqrt(3 * HBAR * GAMMA / (4 * KV ** 3))
np.testing.assert_allclose(u / d_unit, -np.eye(3), atol=1e-14)
np.testing.assert_allclose(v / d_unit, u.conj() @ g / d_unit, atol=1e-14)
print(np.shape(u))
print(u[0])
print(v[0])

outer = np.outer(u[0], v[0]) / abs(u[0, 0]) / abs(v[0, 2])
outer = outer
print(np.round(outer, 2))
print(outer[1])
