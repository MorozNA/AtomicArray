from src.radiative_shift import MarkovianSigmaMatrixForV
from src.radiative_shift import HexagonModel
from src.radiative_shift.dyson_solvers.param import LBAR, GAMMA, OM, KV
import numpy as np

import time
start_time = time.time()

L = 2
DEN = 20

m0 = [-1, 0, 1]
m = [0]

l = L * 2 * np.pi * LBAR
r = 200 / 780 * 2 * np.pi * LBAR
density = DEN * KV ** 3
model = HexagonModel(l, r, density)
model.rotate_model(np.pi/6)

sigma_v = MarkovianSigmaMatrixForV(model)
resolvent = sigma_v.get_resolvent_for_v(OM)
model.add_atom(1.0 * r)

x = np.linspace(1.0, 5.0, 50)
y = np.zeros((len(m), len(x)), dtype='complex_')

for i in range(len(x)):
    radius = r * x[i]
    model.replace_atom(radius)
    s = sigma_v.get_sigma_outside(model, resolvent, m0, m)
    eigs, eigv = np.linalg.eig(s)
    y[:, i] = eigs[:]

y = y / GAMMA
print("--- %s seconds ---" % (time.time() - start_time))

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

# https://stackoverflow.com/questions/11367736/matplotlib-consistent-font-using-latex
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['font.size'] = 12
# mpl.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

L = 2
color = ['violet', 'purple', 'darkblue', 'blue', 'cyan', 'darkgreen', 'lime', 'yellow', 'orange', 'red', 'maroon']
label = ['1', '2', '3', '4']

fig, ax = plt.subplots(1, figsize=(12, 6))

# Figure 1 ____________________________________________

x = np.linspace(1.0, 5.0, len(y[0]))
ax.plot(x, -np.imag(y[0, :]) * 2, linestyle='solid', color=color[3], label='Microscopic approach')

# Labels
ax.set_ylabel(r'$\gamma (r) / \gamma_0$', fontsize=20)

# Asymptote
ax.axhline(y=1, color='black', linewidth=1.0, linestyle='dashed')  # linestyle = (0, (3, 10, 1, 10))

# Axis locators
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.tick_params(axis='both', which='major', direction='in')
ax.tick_params(axis='both', which='minor', direction='in')

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

# Labels
ax.set_xlabel(r'$r / a$', fontsize=20)

# BOTH FIGURES ___________________________________

# Axis limits
ax.set_xlim(1.1, 4.0)

ax.set_ylim(1.0 - 0.1, 1.6)

# Remove weird formatter (+1e2 etc.)
ax.get_yaxis().get_major_formatter().set_useOffset(False)

# Title
plt.suptitle("$^{87}$Ru,  $L = 1 \lambda$, $n_0 \lambdabar^3 = 10$", fontsize=14)

lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

exact_y = [1.6091, 1.4251, 1.2927, 1.1968, 1.1227, 1.0790, 1.0459, 1.0248, 1.0127, 1.0072, 1.0062, 1.0080, 1.0112,
           1.0147, 1.0177, 1.0195, 1.0199, 1.0188, 1.0163, 1.0126, 1.0081, 1.0032, 0.9984, 0.9941, 0.9906, 0.9881,
           0.9869, 0.9870, 0.9882, 0.9904, 0.9933, 0.9967, 1.0001, 1.0033, 1.0060, 1.0080, 1.0091, 1.0092, 1.0085,
           1.0071, 1.0050]
exact_x = np.linspace(1, 5, len(exact_y))

vp_x = exact_x
vpgreen_y = [1.1691, 1.1281, 1.0987, 1.0770, 1.0607, 1.0483, 1.0386, 1.0311, 1.0252, 1.0205, 1.0167, 1.0137, 1.0113,
             1.0093, 1.0077, 1.0063, 1.0053, 1.0044, 1.0036, 1.0030, 1.0025, 1.0021, 1.0018, 1.0015, 1.0012, 1.0010,
             1.0009, 1.0007, 1.0006, 1.0005, 1.0004, 1.0004, 1.0003, 1.0003, 1.0002, 1.0002, 1.0002, 1.0001, 1.0001,
             1.0001, 1.0001]
vpblue_y = [0.3522, 0.2745, 0.2164, 0.1721, 0.1380, 0.1113, 0.0902, 0.0734, 0.0600, 0.0492, 0.0405, 0.0334, 0.0276,
            0.0228, 0.0190, 0.0158, 0.0131, 0.0109, 0.0091, 0.0076, 0.0064, 0.0054, 0.0045, 0.0038, 0.0032, 0.0027,
            0.0022, 0.0019, 0.0016, 0.0013, 0.0011, 0.0010, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0004, 0.0003,
            0.0003, 0.0002]
vpblack_y = [x + y for x, y in zip(vpgreen_y, vpblue_y)]

ax.plot(np.linspace(1, 5, len(exact_y)), exact_y, color='red', label='Fermi\'s golden rule, PRA 97, 023827')
ax.plot(np.linspace(1, 5, len(vpgreen_y)), vpgreen_y, color='green',
        label='Fundamental mode and natural decay diff., PRA 97, 023827')

# Legend
ax.legend(fontsize=16, ncol=1, handleheight=2.4, labelspacing=0.05)

plt.show()
