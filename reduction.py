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
            loss += rate
        elif diag_sp in re_dict[re][1]:
            gain += rate
    if np.max([gain, loss]) == 0:
        print('gain: ', str(gain))
        print('loss: ', str(loss))
        print(re_wM_dict[re])
    return np.max([gain, loss])

def calc_dfj(dat, diag_sp, n_layer, spec_list):
    gain, loss = 0, 0
    gain_old, loss_old = 0, 0
    for re in range(1,nr+1):
        rate = dat['variable']['k'][re][n_layer].astype(float)
        rate_old = np.copy(rate) # no k_time sadly, so rough estimate
        for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
            if sp == 'M': 
                rate *= dat['atm']['n_0'][n_layer]
                rate_old *= dat['atm']['n_0'][n_layer]
            else: 
                rate *= dat['variable']['y_time'][-1, n_layer, spec_list.index(sp)]
                rate_old *= dat['variable']['y_time'][-2, n_layer, spec_list.index(sp)]
        if diag_sp in re_dict[re][0]:
            loss += rate
            loss_old += rate_old
        elif diag_sp in re_dict[re][1]:
            gain += rate
            gain_old += rate_old
    return (gain - loss) - (gain_old - loss_old)

def calc_dfj_dni(dat, spec, diag_sp, n_layer, spec_list):
    dfj = calc_dfj(dat, diag_sp, n_layer, spec_list)
    dni = dat['variable']['y_time'][-1, n_layer, spec_list.index(spec)] - dat['variable']['y_time'][-2, n_layer, spec_list.index(spec)]
    if dni == 0: dni = 1e-10
    return dfj / dni

def calc_sensitivity_param(dat, spec, reduced_list, n_layer):
    B_i = 0
    species = dat['variable']['species']
    for r_sp in reduced_list:
        ni = dat['variable']['y'][n_layer, species.index(spec)]
        gj = calc_gj(dat, r_sp, n_layer, species)
        dfj_dni = calc_dfj_dni(dat, spec, r_sp, n_layer, species)
        B_i += ( (ni/gj)*dfj_dni )**2
    return B_i
        

vul_data = '/scratch/s2555875/output/archean.vul'
with open(vul_data, 'rb') as handle:
    data = pickle.load(handle)
red_spec = ['HCN']
#%%
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
nl = 0
species_list = data['variable']['species']
it = 0 # tracking number of iterations
while again:
    it += 1
    B, B_red = [], []
    for mol in species_list:
        B_new = calc_sensitivity_param(data, mol, red_spec, nl)
        if mol in red_spec:
            B_red.append(B_new)
            B.append(0) # so it will not be added again but index is kept as in species
        else:
            B.append(B_new)
    if  max(B) / min(B_red) > 1e-2:
        red_spec.append(species_list[np.argmax(B)])
    else:
        again = False
        
print('Number of iterations in layer {}: {}'.format(nl, it))
# %%
