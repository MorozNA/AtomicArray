from dataclasses import replace
import numpy as np
from src.radiative_shift import MarkovianSigmaMatrixForV
from src.radiative_shift import CubicComb
from src.radiative_shift.atomspecies import AtomSpecies
from src.radiative_shift.constants import HBAR, C
from src.radiative_shift.tools import find_kd
from draw import draw
from tqdm import tqdm


# 133Cs parameters are F0=4, F=5, J0=1/2, J=3/2, I=7/2
# 87Rb parameters are F0=1, F=0, J0=1/2, J=3/2, I=3/2
F0, F, J0, J, I = 4, 5, 1/2, 3/2, 7/2

N_refr = 3.31

reference_atom = AtomSpecies(F0=F0, F=F, J0=J0, J=J, I=I, lambda_nm=852, gamma=32.815e6)
medium_atom = AtomSpecies(F0=0, F=1, J0=0, J=1, I=0, lambda_nm=852, gamma=reference_atom.gamma)
LBAR, KV = reference_atom.lbar, reference_atom.wavenumber
GAMMA, OM = reference_atom.gamma, reference_atom.omega

a = 2 * np.pi * LBAR / 2 * 0.85147059 # 1.055
height = 0.5 * a
width = 0.5 * a
length_etched = 1.5 * a
length = a * 5
num_etched = int(length / a)
V = (length * height * width + num_etched * (length_etched * height * (a / 2)))
n0 = 50
density = n0 / LBAR ** 3
# num = 2000
# density = num / V
print(int(V * density))
print(int(density * LBAR ** 3))

model = CubicComb(length, a, density, medium_atom, reference_atom)
print(len(model.x))
print(len(model.x) / V * LBAR ** 3)

# Use the comb volume already calculated above for the actual density.
detuning = find_kd(N_refr, len(model.x) / V * LBAR ** 3)
medium_omega = OM - detuning * medium_atom.gamma
model.medium_atom = replace(medium_atom, lambda_nm=2 * np.pi * C / medium_omega * 1e7)
model._refresh_properties()
# The solvers still use the old distance-method name.
model.calculate_distances_to_signal_atom = model.calculate_distances_to_reference_atom

# Quantize along the original x axis.
model.rotate_about_y(-np.pi/2)
sigma = MarkovianSigmaMatrixForV(model)
resolvent = sigma.get_resolvent_for_v(OM)
model.rotate_about_y(np.pi/2)

model.set_reference_position_cylindrical(reference_atom)

x = np.linspace(0.01, 5.0, 200)
y = np.zeros((len(reference_atom.m), len(x)), dtype=complex)

for i in range(len(x)):
    model.set_reference_position_cylindrical(reference_atom, 200e-7 * x[i])

    model.rotate_about_y(-np.pi/2)
    s = sigma.get_sigma_outside(model, resolvent)
    model.rotate_about_y(np.pi/2)

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

    if i == 10:
        eigv_saved = np.round(eigv, 3)
        eigs_saved = y[:, i]
        print('Saving eigenvectors...')
        print(i)
        print(x[i])
        print(eigs_saved / HBAR / GAMMA)
        for index in range(len(reference_atom.m)):
            print(f'eigenvalue {index}: ', np.round(eigv_saved[:, index], 1))

y = y / HBAR / GAMMA
np.savetxt('./data/fig5.txt', y, fmt='%f')
np.savetxt('./data/fig5_eigv40.txt', eigv_saved, fmt='%f')
draw(x, y, 10)
