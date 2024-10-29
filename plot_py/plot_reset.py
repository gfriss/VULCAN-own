import matplotlib.pyplot as plt
from cycler import cycler

def reset_plt(ticksize, fontsize, fxsize, fysize):
    plt.style.use('seaborn-white') # was seaborn-v0_8-white
    plt.rcParams['xtick.labelsize'] = ticksize
    plt.rcParams['ytick.labelsize'] = ticksize
    plt.rcParams['xtick.major.width'] = 1.25
    plt.rcParams['ytick.major.width'] = 1.25
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['xtick.major.size'] = 4
    plt.rcParams['ytick.major.size'] = 4
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['ytick.minor.size'] = 4
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['mathtext.fontset'] = 'dejavuserif' #'stix'
    plt.rcParams['font.family'] = 'serif' #'STIXGeneral'
    plt.rcParams['legend.facecolor'] = 'white'
    plt.rcParams['legend.fontsize'] = ticksize
    plt.rcParams['axes.formatter.limits'] = (-2,4)
    plt.rcParams['axes.linewidth'] = 1.75
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.6
    plt.rcParams['figure.figsize'] = (fxsize,fysize)
    plt.rcParams['figure.autolayout'] = True
    cm = plt.get_cmap('nipy_spectral')
    n_colours = 16
    lst = ['-', '--', '-.', ':']
    plt.rcParams['axes.prop_cycle'] = cycler(linestyle = n_colours//len(lst) * lst) + cycler(color = [cm(1.*i/n_colours) for i in range(n_colours)])