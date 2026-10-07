import numpy as np
from two_atoms_calc_tcs import calc_tsc, LBAR


tcs1 = calc_tsc([0, 0, 0], [0, 0, 0.5 * LBAR])
tcs2 = calc_tsc([0, 0, 0], [0, 0, 1.0 * LBAR])
tcs3 = calc_tsc([0, 0, 0], [0, 0, 10.0 * LBAR])


from matplotlib import pyplot as plt
x = np.linspace(-15, 15, 1000)
plt.plot(x, tcs1, '-', color='tab:red', label=r'$r / \bar\lambda = 0.5$')
plt.plot(x, tcs2, color='tab:blue', label=r'$r / \bar\lambda = 1$')
plt.plot(x, tcs3, color='tab:grey', label=r'$r / \bar\lambda = 10$')
plt.xlabel(r'$\Delta / \gamma$')
plt.ylabel(r'$\sigma_{\mathrm{tot}} / (\lambda / 2\pi)^2$')
plt.legend()
plt.show()
