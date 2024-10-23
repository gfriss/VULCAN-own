#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pickle
import plot_reset as pr
import os
wd = os.getcwd()
os.chdir('../')
from chem_funs import nr, re_wM_dict, re_dict
os.chdir(wd)
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6)

out_folder = '/scratch/s2555875/output/'
plot_folder = '/scratch/s2555875/plot/'

vul_data = '/scratch/s2555875/output/archean.vul'
with open(vul_data, 'rb') as handle:
    data = pickle.load(handle)
species = data['variable']['species']
#%%
def get_total_reaction_rate_evol(dat, diag_sp):
    t = dat['variable']['t_time']
    total_re_list = np.zeros((len(t), 120))
    for re in range(1,nr+1):
        if diag_sp in re_wM_dict[re][0] or diag_sp in re_wM_dict[re][1]:
            rate = dat['variable']['k'][re].astype(float) # rate coefficient does not change with time (it's calculated from the Arrhenius formula)
            for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
                if sp == 'M': rate *= dat['atm']['n_0']
                else: rate = dat['variable']['y_time'][:, :, species.index(sp)] * rate
            if diag_sp in re_dict[re][0]:
                total_re_list -= np.array(rate)
            elif diag_sp in re_dict[re][1]:
                total_re_list += np.array(rate)
    return total_re_list    

def plot_reaction_rate_evol_layer(t, k, layer, figname = None):
    fig, ax =plt.subplots()
    ax.plot(t, k[:, layer])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'$k_{tot}$ [cm$^3$s$^{-1}$]')
    ax.set_xscale('log')
    ax.set_yscale('symlog')
    
    if figname != None:
        fig.savefig(plot_folder + figname)
        
def log_tick_formatter(val, pos = 0):
    return f"$10^{{{int(val)}}}$"

def symlog_tick_formatter(val, pos = 0):
    if val > 0:
        return f"$10^{{{int(val)}}}$"
    else:
        return f"$-10^{{{int(abs(val))}}}$"
  
def plot_reaction_rate_evol_all(t, k, figname = None):
    layers = np.arange(0, 120, 1)
    L, T = np.meshgrid(layers, np.log10(t))
    k[k > 0] = np.log10(k[k > 0])
    k[k < 0] = -np.log10(abs(k[k < 0]))
    fig, ax = plt.subplots(subplot_kw = {'projection': '3d'})
    ax.plot_wireframe(L, T, k, rstride=10, cstride=10, color = 'k', alpha = 0.6)
    #for l in layers:
    #    ax.plot(np.ones_like(t)*l, np.log10(t), k[:, l], c = 'k', alpha = 0.6)
    ax.set_xlabel('# layers')
    ax.invert_xaxis()
    ax.set_ylabel('Time [s]')
    ax.set_zlabel(r'$k_{tot}$ [cm$^3$s$^{-1}$]')
    ax.view_init(elev=20., azim=-55, roll = 0)
    ax.zaxis.set_major_formatter(mticker.FuncFormatter(symlog_tick_formatter))
    ax.zaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    if figname != None:
        fig.savefig(plot_folder + figname)
#%%
spec = 'H2CO'
k_layers = get_total_reaction_rate_evol(data, spec)
#%%
nl = 0
plot_reaction_rate_evol_layer(data['variable']['t_time'], k_layers, nl)
# %%
plot_reaction_rate_evol_all(data['variable']['t_time'], k_layers)
# %%
