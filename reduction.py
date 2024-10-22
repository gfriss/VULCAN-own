#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
from chem_funs import nr, re_wM_dict, re_dict
import plot_py.plot_reset as pr
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6)

def calc_gj(dat, diag_sp, n_layer, spec_list):
    gain, loss = 0, 0
    for re in range(1,nr+1):
        rate = dat['variable']['k'][re][n_layer].astype(float)
        for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
            if sp == 'M': rate *= dat['atm']['n_0'][n_layer]
            else: rate *= dat['variable']['y'][n_layer, spec_list.index(sp)]
        if diag_sp in re_dict[re][0]:
            loss += abs(rate)
        elif diag_sp in re_dict[re][1]:
            gain += abs(rate)
    #if np.max([gain, loss]) == 0:
    #    print('gain: ', str(gain))
    #    print('loss: ', str(loss))
    #    print(re_wM_dict[re])
    return np.max([gain, loss])

def calc_dfj(dat, spec, diag_sp, n_layer, spec_list, dni):
    gain, loss = 0, 0
    gain_delta, loss_delta = 0, 0
    for re in range(1,nr+1):
        rate = dat['variable']['k'][re][n_layer].astype(float)
        rate_delta = np.copy(rate)
        for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
            if sp == 'M': 
                rate *= dat['atm']['n_0'][n_layer]
                rate_delta *= dat['atm']['n_0'][n_layer]
            else: 
                rate *= dat['variable']['y'][n_layer, spec_list.index(sp)]
                if sp == spec:
                    rate_delta *= dat['variable']['y'][n_layer, spec_list.index(sp)] + dni # only veried for species i
                else:
                    rate_delta *= dat['variable']['y'][n_layer, spec_list.index(sp)]
        if diag_sp in re_dict[re][0]:
            loss_delta += abs(rate_delta)
            loss_delta += abs(rate_delta)
        elif diag_sp in re_dict[re][1]:
            gain_delta += abs(rate_delta)
            gain_delta += abs(rate_delta)
    return (gain_delta - loss_delta) - (gain - loss)

def calc_dfj_dni(dat, spec, diag_sp, n_layer, spec_list):
    #dfj = calc_dfj(dat, diag_sp, n_layer, spec_list)
    dni = dat['variable']['y_time'][-1, n_layer, spec_list.index(spec)] - dat['variable']['y_time'][-2, n_layer, spec_list.index(spec)]
    if dni == 0: # no dependency
        return 0
    else:
        dfj = calc_dfj(dat, spec, diag_sp, n_layer, spec_list, dni)
        return dfj / dni

def calc_sensitivity_param(dat, spec, reduced_list, n_layer):
    B_i = 0
    species = dat['variable']['species']
    for r_sp in reduced_list:
        ni = dat['variable']['y'][n_layer, species.index(spec)]
        gj = calc_gj(dat, r_sp, n_layer, species)
        if gj != 0: # if this is zero it means that the species is not present at that layer
            dfj_dni = calc_dfj_dni(dat, spec, r_sp, n_layer, species)
            B_i += ( (ni/gj)*dfj_dni )**2
    return B_i
        
def get_rates_and_sums(dat, n_layer, spec_list):
    sum_rate_per_species = np.zeros(len(spec_list))
    rate_list = []
    for re in range(1,nr+1):
        rate = dat['variable']['k'][re][n_layer].astype(float)
        sp_idx_per_reaction = []
        for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
            if sp == 'M': rate *= dat['atm']['n_0'][n_layer]
            else: 
                rate *= dat['variable']['y'][n_layer, spec_list.index(sp)]
                sp_idx_per_reaction.append(spec_list.index(sp))
        for sp in re_wM_dict[re][1]: # [0]: reactants; [1]: prodcuts
            if sp != 'M':
                sp_idx_per_reaction.append(spec_list.index(sp))
        rate_list.append(rate)
        for idx in sp_idx_per_reaction:
            sum_rate_per_species[idx] += rate # will it be repeated?
    return rate_list, sum_rate_per_species
            
def calc_wr(wr, ws, rate_list, sum_rate_per_species, iter_specs, reduced_list):
    spec_to_add = []
    for i,sp in enumerate(iter_specs):
        for re in range(1,nr+1):
            if (sp in re_wM_dict[re][0] or sp in re_wM_dict[re][1]) and sp != 'M':
                wr_new = (rate_list[re-1] / sum_rate_per_species[i]) * ws[i]
                wr[re-1] = np.max([wr[re-1], wr_new])
                spec_to_add += [s for s in (list(re_wM_dict[re][0])+list(re_wM_dict[re][1])) if s not in reduced_list and s != 'M']
    spec_to_add = np.unique(spec_to_add) # making sure that no duplicates are present
    for i,sp in enumerate(spec_to_add):
        wr_i_max = 0
        for re in range(1,nr+1):
            if sp in re_wM_dict[re][0] or sp in re_wM_dict[re][1]:
                wr_i_max = np.max([wr_i_max, wr[re-1]]) # finding maximum reaction weith coupled to current species
        ws[i] = np.max([ws[i], wr_i_max])
    return wr, ws, list(spec_to_add)

vul_data = '/scratch/s2555875/output/archean.vul'
with open(vul_data, 'rb') as handle:
    data = pickle.load(handle)
#%%
#red_spec = ['HCN']
#B = []
#for mol in data['variable']['species']:
#    if mol not in red_spec:
#        B.append(calc_sensitivity_param(data, mol, red_spec, 0))
#    else:
#        B.append(0) # set to zero so it would not be given to the list again
        
#%% 
#B.sort(reverse = True)
#plt.plot(B, linestyle = '', marker = '.')
#plt.yscale('log')
#%%
again = True
nl = 100
species_list = data['variable']['species']
it = 0 # tracking number of iterations
red_spec = ['HCN', 'HCN_rain', 'H2O_rain']
B_test = []
while again:
    it += 1
    B, B_red, B_out = [], [], []
    for mol in species_list:
        B_new = calc_sensitivity_param(data, mol, red_spec, nl)
        B.append(B_new)
        if mol in red_spec:
            B_red.append(B_new)
            B_out.append(0) # to keep its indexing and not change the max value and the total B
        else:
            B_out.append(B_new)
    if  max(B_out) / min(B_red) > 1e-2:
        red_spec.append(species_list[np.argmax(B_out)])
    else:
        again = False
        
print('Number of iterations in layer {}: {}'.format(nl, it))
B_sort_idx = np.argsort(B)
B.sort(reverse = True)
plt.plot(B, linestyle = '', marker = '.', color = 'r')
plt.yscale('log')
# %%
nl = 113
species_list = data['variable']['species']
rates, sum_rates = get_rates_and_sums(data, nl, species_list)
it = 0 # tracking number of iterations
red_spec = ['HCN', 'HCN_rain', 'H2O_rain']
again = True
w_r = np.zeros_like(rates)
w_s = np.zeros(len(species_list))
for rs in red_spec: w_s[species_list.index(rs)] = 1.
while again:
    it += 1
    if it == 1: 
        added_spec = np.copy(red_spec)
    w_r, w_s, added_spec = calc_wr(w_r, w_s, rates, sum_rates, added_spec, red_spec)
    red_spec += added_spec
    if len(red_spec) == len(species_list):
        again = False
        
w_r_sort_idx = np.argsort(w_r)
w_r = list(w_r)
w_r.sort(reverse = True)
plt.plot(w_r, linestyle = '', marker = '.', color = 'r')
plt.xlabel('# reaction')
plt.ylabel(r'w_r')
plt.yscale('log')
plt.ylim((1e-18, None))
plt.savefig('/scratch/s2555875/plot/archean_reduction_{}.pdf'.format(nl))
with open('/scratch/s2555875/plot/archean_reduced_network.txt', 'a') as f:
    f.write('Layer {}:\n'.format(nl))
    for i,re in enumerate(w_r_sort_idx[np.array(w_r) > 1e-18]):
        f.write(str(re_wM_dict[re+1]) + ' - w_r = {}\n'.format(w_r[i])) # maybe need to change and chek if re%2 == 0...
    f.write('\n') # just to space out things
# %%
