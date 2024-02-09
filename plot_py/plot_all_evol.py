import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pickle

vul_file = sys.argv[1]

with open(vul_file, 'rb') as handle:
    data = pickle.load(handle)

vulcan_spec = data['variable']['species']

def plot_time_evo(yt, tt, n, mol, xscale = 'linear', yscale = 'log', savefig = False):
    alpha = np.linspace(0.1, 1, n)

    fig, ax = plt.subplots()
    for layer in range(n):
        if layer in [0, n//4, n//2, 3*n//4, n-1]:
            ax.plot(tt, yt[:, layer, vulcan_spec.index(mol)], label = 'layer = {}'.format(layer), alpha = alpha[layer], color = 'red')
        else:
            ax.plot(tt, yt[:, layer, vulcan_spec.index(mol)], alpha = alpha[layer], color = 'red')

    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'n [cm$^{-3}$]')
    ax.set_title(mol)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.legend()
    ax.set_ylim((np.min(yt[-1, :, vulcan_spec.index(mol)]), None))
    if savefig:
        fig.savefig('../plot/all_evol/{}.png'.format(mol))
    plt.close(fig)
# %%
def plot_evo_layer(yt, tt, n, mol, xscale = 'linear', yscale = 'log', ylim = None):
    fig, ax = plt.subplots()

    ax.plot(tt, yt[:, n, vulcan_spec.index(mol)], color = 'red', label = mol)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'n [cm$^{-3}$]')
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.legend()
    if ylim != None:
        ax.set_ylim(ylim)

for species in vulcan_spec:
    plot_time_evo(data['variable']['y_time'], data['variable']['t_time'], 120, species, savefig = True)