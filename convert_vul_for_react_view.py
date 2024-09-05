import numpy as np
import pickle

# ------------creating the species_comp.out file------------
allcompose = np.genfromtxt('thermo/all_compose.txt', dtype = str)
phase = 'Gas'
surf_cov = '0'
array_insert = np.array([phase, surf_cov])

new_str = '{: <15}{: ^7}{: ^9}'.format('Species_name', 'Phase', 'Surf_conv')
for element in allcompose[0, 1:-1]: # excluding species and mass
    new_str += '\t{}'.format(element)
new_str += '\n'

for line in allcompose[1:, :-1]:
    new_str += '{: <15}{: ^7}{: ^9}\t'.format(line[0], array_insert[0], array_insert[1])
    new_str += '\t'.join(line[1:]) + '\n'

with open('/scratch/s2555875/species_comp.out', 'w') as f:
    f.write(new_str)

# ------------creating the reaction_rate.out file------------
from chem_funs import nr, re_wM_dict

new_str_reaction = 'Fwd_Rate\tRev_Rate\tNet_Rate\tPEI\tReaction_String\n'

with open('/scratch/s2555875/output/archean.vul', 'rb') as handle: # need to generalise
    data = pickle.load(handle)
species = data['variable']['species']
n = 42 # need to generalise
for re in range(1, nr, 2):
    fwd_rate, rev_rate = 0., 0.
    rate = data['variable']['k'][re][n].astype(float)
    for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
        if sp == 'M': rate *= data['atm']['n_0'][n]
        #else: rate *= data['variable']['y'][:,species.index(sp)] 
        elif sp != 'S8_l_s': rate *= data['variable']['y'][:,species.index(sp)][n] 
        fwd_rate = rate
    for sp in re_wM_dict[re+1][0]: # [0]: reactants; [1]: prodcuts
        if sp == 'M': rate *= data['atm']['n_0'][n]
        #else: rate *= data['variable']['y'][:,species.index(sp)] 
        elif sp != 'S8_l_s': rate *= data['variable']['y'][:,species.index(sp)][n] 
        rev_rate = rate
    net_rate = fwd_rate - rev_rate
    if fwd_rate == 0 and rev_rate == 0:
        pei = 0.
    else:
        pei = fwd_rate / (fwd_rate + rev_rate)
    reaction_str = ' + '.join(re_wM_dict[re][0]) + ' <=> ' + ' + '.join(re_wM_dict[re][1])
    new_str_reaction += '{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{}\n'.format(fwd_rate, rev_rate, net_rate, pei, reaction_str)

with open('/scratch/s2555875/reaction_rates.out', 'w') as f:
    f.write(new_str_reaction)