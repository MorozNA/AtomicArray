from src.radiative_shift import VMediumSelfEnergyMatrix
from src.radiative_shift import HexagonSphere
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import HBAR
import numpy as np

atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=780, gamma=38.11e6)
# Keep the dense eigensystems small enough for an interactive test.
RADIUS = atom.lbar
density = 10 / atom.lbar ** 3

x = np.linspace(0.1, 1.5, 10)
y = []
for i in x:
    radius = RADIUS * i
    model = HexagonSphere(radius, density, atom, atom)
    natoms = len(model.x)
    dens = 3 * natoms / (4 * np.pi * radius**3)
    print(f"Density = {dens * atom.lbar ** 3:.3f} n lambda_bar^3; atoms = {natoms}")

    sigma = VMediumSelfEnergyMatrix(model)
    eigs = np.linalg.eigvals(sigma.sigma / (HBAR * atom.gamma))
    eigs = np.sort(np.imag(eigs))
    # Vacuum decay is included in the resolvent, not in sigma.sigma.
    decay_rates = 1 - 2 * eigs
    assert np.isfinite(eigs).all()
    assert np.min(decay_rates) >= -1e-10
    np.testing.assert_allclose(np.sum(decay_rates), 3 * natoms, rtol=1e-12)
    print("Minimum decay / gamma = {:.5g}".format(decay_rates[-1]))
    print("Maximum decay / gamma = {:.5g}".format(decay_rates[0]))
    y.append(-(radius / atom.lbar)**2 * eigs[0] / natoms)  # see PRA 93, 043830

from matplotlib import pyplot as plt
plt.plot(x,y)
plt.xlabel(r'$R / \bar\lambda$')
plt.ylabel(r'$-(R / \bar\lambda)^2\,\mathrm{Im}\,\Sigma_{\min} / (N\hbar\gamma)$')
plt.show()
