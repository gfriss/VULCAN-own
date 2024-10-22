#%%
import numpy as np
import matplotlib.pyplot as plt
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
def get_total_reaction_rate_layer(dat, diag_sp, layer):
    t = dat['variable']['t_time']
    total_re_list = np.zeros_like(t)
    for re in range(1,nr+1):
        if diag_sp in re_wM_dict[re][0] or diag_sp in re_wM_dict[re][1]:
            rate = dat['variable']['k_time'][re, layer].astype(float)
            for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
                if sp == 'M': rate *= dat['atm']['n_0']
                else: rate *= dat['variable']['y_time'][:, layer, species.index(sp)]
            if diag_sp in re_dict[re][0]:
                total_re_list -= np.array(rate)
            elif diag_sp in re_dict[re][1]:
                total_re_list += np.array(rate)
    return total_re_list    
#%%
nl = 0
spec = 'HCN'
k_layer = get_total_reaction_rate_layer(data, spec, nl)
fig, ax =plt.subplots()
ax.plot(data['variable']['t_time'], k_layer)
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'$k_{tot}$ [cm$^3$s$^{-1}$]')