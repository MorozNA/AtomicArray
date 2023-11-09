import numpy as np
from two_atoms_calc_tcs import calc_tsc
from src.radiative_shift.dyson_solvers.param import LBAR


tcs1 = calc_tsc([0, 0, 0], [0, 0, 0.5 * LBAR])
tcs2 = calc_tsc([0, 0, 0], [0, 0, 1.0 * LBAR])
tcs3 = calc_tsc([0, 0, 0], [0, 0, 10.0 * LBAR])


from matplotlib import pyplot as plt
x = np.linspace(-15, 15, 1000)
plt.plot(x, tcs1, '-', color='tab:red')
plt.plot(x, tcs2, color='tab:blue')
plt.plot(x, tcs3, color='tab:grey')
plt.xlabel('$\Delta / \gamma$')
plt.ylabel('$\sigma_0$ in units of $\lambda / 2 \pi^3$')
plt.show()
