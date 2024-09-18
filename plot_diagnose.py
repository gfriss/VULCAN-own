# Diagnostic tool for output largest rates
#%%
import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle

import chem_funs, op
from chem_funs import nr, re_wM_dict, re_dict

vul_data = '/scratch/s2555875/output/archean.vul'

# setting the numerical solver to the desinated one in vulcan_cfg
solver_str = vulcan_cfg.ode_solver
solver = getattr(op, solver_str)()

# the number of fastest reactions to print out
top_num = 205

# the species involved
diag_sp = 'HCN'

# the region to look at
# the bottom layer 
n_bot = 0
# the top layer
n_top = 120

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
species = data['variable']['species'] 
pressure = data['atm']['pco']/1e6 # bar
#%%
rate_list = []
max_re_list = [] 
total_re_list = np.zeros_like(pressure) 
for re in range(1,nr+1):
    #if diag_sp in re_wM_dict[re][0] or diag_sp in re_wM_dict[re][1]:
    rate = data['variable']['k'][re][n_bot:n_top].astype(float)
    for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
        if sp == 'M': rate *= data['atm']['n_0'][n_bot:n_top]
        else: rate *= data['variable']['y'][:,species.index(sp)][n_bot:n_top] 
    rate_list.append(rate)
    max_re_list.append(np.amax(rate))
    if diag_sp in re_dict[re][0]:
        total_re_list -= np.array(rate)
    elif diag_sp in re_dict[re][1]:
        total_re_list += np.array(rate)

# 1+ is to shift the index to match starting with R1  
re_sort_indx = 1 + np.argsort(max_re_list)[::-1]
#rate_sort = np.sort(max_re_list)[::-1]
top_re = re_sort_indx[0:top_num]
#%%
fig, ax = plt.subplots()
for re in top_re:
    if diag_sp in re_dict[re][0] or diag_sp in re_dict[re][1]:
        if re % 2 == 1: 
            ax.plot(rate_list[re-1], pressure, label = str(data['variable']['Rf'][re]))
        else: 
            ax.plot(rate_list[re-1], pressure, label = 'rev '+str(data['variable']['Rf'][re-1]))
ax.plot(total_re_list, pressure, 'k--', label = 'total rate')

ax.invert_yaxis()
ax.legend(bbox_to_anchor = (1.05, 0.85))
ax.set_xlabel('Reaction rate [cm3s-1]')
ax.set_ylabel('Pressure [bar]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-15,None)
# %%
fig, ax = plt.subplots()
for re in top_re:
    if diag_sp in re_dict[re][0]:
        if re % 2 == 1: 
            ax.plot(rate_list[re-1], pressure, label = str(data['variable']['Rf'][re]), linestyle = '-')
        else: 
            ax.plot(rate_list[re-1], pressure, label = 'rev '+str(data['variable']['Rf'][re-1]), linestyle = '-')
    elif diag_sp in re_dict[re][1]:
        if re % 2 == 1: 
            ax.plot(rate_list[re-1], pressure, label = str(data['variable']['Rf'][re]), linestyle = '--')
        else: 
            ax.plot(rate_list[re-1], pressure, label = 'rev '+str(data['variable']['Rf'][re-1]), linestyle = '--')
ax.plot(total_re_list, pressure, 'k:', label = 'total rate')

ax.invert_yaxis()
ax.legend(bbox_to_anchor = (1.05, 0.85))
ax.set_xlabel('Reaction rate [cm3s-1]')
ax.set_ylabel('Pressure [bar]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-15,None)
# %%
