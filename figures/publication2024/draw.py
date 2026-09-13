from src.radiative_shift import EmptyModel
import numpy as np
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


def draw(x, y, point=0, title='', ind=None, cond=0, monte_x=None, monte_y=None):
    if ind is None:
        ind = [x for x in range(y.shape[0])]

    if cond == 0:
        cond = [True] * y.shape[0]
    elif cond == 1:
        cond = [ind[i] > 0 or ind[i] == 0 for i in range(y.shape[0])]
    elif cond == 2:
        cond = [ind[i] < 0 or ind[i] == 0 for i in range(y.shape[0])]
        print(cond)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 6), gridspec_kw={'hspace': 0})
    for i in range(y.shape[0]):
        if cond[i]:
            ax[0].plot(x[point:], - 2 * np.imag(y[i, point:]), color=color[i], label=str(ind[i]))
            ax[1].plot(x[point:], np.real(y[i, point:]), color=color[i])
    ax[0].set_ylabel(r'$\Gamma (x) / \Gamma(\infty)$')
    ax[1].set_ylabel(r'$\Delta (x) / \Gamma(\infty)$')
    ax[1].set_xlabel(r'$x / x_0$')

    if monte_x is None:
        pass
    else:
        ax[1].plot(monte_x, monte_y, linestyle='dashed', color='black')

    # Asymptote
    ax[0].axhline(y=1, color='black', linewidth=1.0, linestyle='dashed')  # linestyle = (0, (3, 10, 1, 10))

    # Axis locators
    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax[0].tick_params(axis='both', which='major', direction='in')
    ax[0].tick_params(axis='both', which='minor', direction='in')

    # Hide the right and top spines
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax[0].yaxis.set_ticks_position('left')
    ax[0].xaxis.set_ticks_position('bottom')

    ax[0].set_xlim(x[point], 2.99)
    ax[1].set_xlim(x[point], 3.0)
    # ax[0].set_ylim(0.9, np.amax(np.imag(y)[:, point:]))
    ax[1].set_ylim(-0.15, 0.15)
    ax[0].set_title(title, fontsize=16)
    ax[0].legend(ncol=2)

    plt.show()

    # TODO: set ylims
