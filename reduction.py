#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
from chem_funs import nr, re_wM_dict, re_dict

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
    sum_rate_per_species = np.zeros_like(spec_list)
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
            
def calc_wr(rate_list, sum_rate_per_species, spec_list, reduced_list, first_iter):
    wr = np.zeros_like(rate_list)
    ws = np.zeros_like(spec_list)
    for i,sp in reduced_list:
        if first_iter:
            ws[i] = 1
        for re in range(1,nr+1):
            wr_new = (rate_list[re] / sum_rate_per_species[i]) * ws[i]
            wr[re] = np.max([wr[re], wr_new])     
            ws[i] = np.max([ws[i], rate_list[re]])   # not ready, to do: keep original species at 1 and only update the new ones
    return wr, ws

vul_data = '/scratch/s2555875/output/archean.vul'
with open(vul_data, 'rb') as handle:
    data = pickle.load(handle)
#%%
red_spec = ['HCN']
B = []
for mol in data['variable']['species']:
    if mol not in red_spec:
        B.append(calc_sensitivity_param(data, mol, red_spec, 0))
    else:
        B.append(0) # set to zero so it would not be given to the list again
        
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
