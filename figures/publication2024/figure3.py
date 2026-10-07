from dataclasses import replace
import numpy as np
from src.radiative_shift import VMediumSelfEnergyMatrix
from src.radiative_shift import HexagonModel
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import HBAR, C
from src.radiative_shift.tools import find_reference_detuning
from draw import draw
from tqdm import tqdm

# 133Cs parameters are F0=4, F=5, J0=1/2, J=3/2, I=7/2
# 87Rb parameters are F0=1, F=0, J0=1/2, J=3/2, I=3/2
F0, F, J0, J, I = 4, 5, 1/2, 3/2, 7/2

L = 2
DEN = 20
N_refr = 1.45

reference_atom = AtomSpecies(F0=F0, F=F, J0=J0, J=J, I=I, lambda_nm=852, gamma=32.815e6)
medium_atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=852, gamma=reference_atom.gamma)
LBAR, KV = reference_atom.lbar, reference_atom.wavenumber
GAMMA, OM = reference_atom.gamma, reference_atom.omega

l = L * 2 * np.pi * LBAR
r = 200 / 852 * 2 * np.pi * LBAR
density = DEN * KV ** 3

model = HexagonModel(l, r, density, medium_atom, reference_atom)
model.set_reference_position_cylindrical(reference_atom, 1.0 * r)
detuning = find_reference_detuning(N_refr, model.properties.density)
medium_omega = OM - detuning * medium_atom.gamma
model.medium_atom = replace(medium_atom, lambda_nm=2 * np.pi * C / medium_omega * 1e7)
model._refresh_properties()

# Quantize along the original x axis.
model.rotate_about_y(-np.pi/2)
sigma = VMediumSelfEnergyMatrix(model)
resolvent = sigma.get_medium_resolvent(OM)


x = np.linspace(1.0, 5.0, 50)
y = np.zeros((len(reference_atom.m), len(x)), dtype=complex)

for i in tqdm(range(len(x))):
    model.rotate_about_y(np.pi/2)
    model.set_reference_position_cylindrical(reference_atom, r * x[i])
    model.rotate_about_y(-np.pi/2)

    s = sigma.get_reference_self_energy(model, resolvent)
    eigs_temp, eigv_temp = np.linalg.eig(s)

    eigv2 = np.zeros([len(reference_atom.m), len(reference_atom.m)], dtype=complex)
    prod = np.zeros((len(reference_atom.m), 1))

    if i > 0:
        for j in range(len(reference_atom.m)):
            for k in range(len(reference_atom.m)):
                prod[k] = abs(np.vdot(eigv[:, j], eigv_temp[:, k]))
            eigv2[:, j] = eigv_temp[:, np.argmax(prod)]
            y[j, i] = (eigs_temp[np.argmax(prod)])
            eigv_temp[:, np.argmax(prod)] = np.zeros(np.shape(eigv_temp[:, 0]))
        eigv = eigv2
    else:
        for j in range(len(reference_atom.m)):
            y[j, i] = eigs_temp[j]
        eigv = np.copy(eigv_temp)
    if i==45:
        for j in range(11):
            print("eigv", j, ":", np.round(eigv[:, j], 2))

y = y / HBAR / GAMMA
np.savetxt('./data/fig3.txt', y, fmt='%f')
draw(x, y, point=1)
