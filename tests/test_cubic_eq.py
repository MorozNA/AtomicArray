from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from src.radiative_shift import find_reference_detuning
from src.radiative_shift.tools.cubic_equation import solve_cubic


n_refr = 1.45
# n_refr = 4.5

F = 1
F0 = 0
J = 1 / 2
J0 = 3 / 2
# n0 = 133.76963635370498
n0 = 10

alpha = 1  # (2 * F0 + 1) / 3
ro = n0 * (2 * F + 1) / 3 / (2 * F0 + 1)

del_omega = np.linspace(-100 * ro, 1 * ro, 10000)  # del to gamma
eps = []

for k in range(len(del_omega)):
    delt = del_omega[k]

    a = 1j / 2
    b = alpha * np.pi * ro + delt
    c = - 1j / 2
    d = alpha * 2 * np.pi * ro - delt

    x = solve_cubic(a, b, c, d)
    eps.append(x[0] ** 2)

y = n_refr ** 2 * np.ones(np.size(del_omega))

from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=[8, 5])

ax.plot(del_omega, np.real(eps), 'r', label=r'solution for $\epsilon(\omega)$')
ax.plot(del_omega, abs(np.imag(eps)), 'r--', label='imaginary part')
ax.set_xlabel(r'$(\omega - \omega_d)/\gamma$', fontsize=18)
ax.set_ylabel(r'$\epsilon (\omega)$', fontsize=18)
ax.plot(del_omega, y, 'b--', label=r'$\epsilon_{Si} = 1.45^2$')
ax.legend(fontsize=14)
ax.set_title(r'$n_0$ = {first:1.5} $\bar\lambda^3 (2F + 1) / [3(2F0 + 1)]$'.format(
    first=ro / (2 * F + 1) * 3 * (2 * F0 + 1)), fontsize=18)  # (2F + 1) / [3(2F0 + 1)]

# Limits

n = n_refr
n1 = np.sqrt(n**2 - 0.05)
n2 = np.sqrt(n**2 + 0.05)

pr_d = np.pi * ro * (2 + n ** 2) / (1 - n ** 2)
pr_d1 = np.pi * (2 + n1 ** 2) / (1 - n1 ** 2)
pr_d2 = np.pi * (2 + n2 ** 2) / (1 - n2 ** 2)

k_del = find_reference_detuning(n_refr, n0)
roots = solve_cubic(1j / 2, alpha * np.pi * ro + k_del,
                    -1j / 2, alpha * 2 * np.pi * ro - k_del)
np.testing.assert_allclose(np.real(roots[0] ** 2), n_refr ** 2, atol=0.01)
print(f"Detuning = {k_del:.6f} gamma; permittivity = {roots[0] ** 2}")
ax.axvline(k_del)
ax.set_xlim(k_del - 1, k_del + 1)
ax.set_ylim(n**2 - 0.1, n**2 + 0.1)

plt.show()
