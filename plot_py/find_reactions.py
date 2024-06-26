#%%
import numpy as np
import pickle
from collections import defaultdict
import matplotlib.pyplot as plt
import plot_reset as ps

def calc_rate(dat, spec_list, re_id, n):
    ''' Calculates the reaction rate for a reaction that is given by
        its reaction id (re_id) using the reaction coefficient (k) and
        number density of reagent species (except third body M).'''
    rate = dat['variable']['k'][re_id][n]
    for sp in spec_list:
        if sp != 'M':
            rate *= dat['variable']['y'][n, dat['variable']['species'].index(sp)]
    return rate

def get_species(eq_side):
    ''' Takes products or reagents string (A + B form) and turns it into a
        list of strings of the species.'''
    side_split = eq_side.split('+')
    if len(side_split) == 1: # stripping them from white spaces
        side_split = np.array([side_split[0].strip()]) # with array length 1 the other method fails so doing it separately
    else:
        side_split = np.array([r.strip() for r in side_split])
    return side_split

def max_re_one_step_layer(dat, mol, n, prod):
    ''' This functions finds the most important production/destruction reaction for a given molecule
        in a given layer and returns the reaction rate, reagents/products and the reaction itself.
    '''
    reaction = ''
    mol_returned = []
    k_rea = 0
    for k,v in dat['variable']['Rf'].items():
        reagents_spec, products_spec = [], []
        re_rate = 0
        reagents_products = v.split('->') # separating production and destruction
        reagents_spec = get_species(reagents_products[0])
        products_spec = get_species(reagents_products[1])
        mol_current = np.array([])
        if prod == False and mol in reagents_spec: # destruction so it is on the reagents side
            re_rate = calc_rate(dat, reagents_spec, k, n)
            mol_current = products_spec
        elif prod == True and mol in products_spec: # production side
            re_rate = calc_rate(dat, products_spec, k, n)
            mol_current = reagents_spec
        if re_rate > k_rea:
            k_rea = re_rate
            reaction = v
            mol_returned = list(mol_current)
    return k_rea, mol_returned, reaction

def walker(dat, mol, n):
    ''' Starting from a given molecule (mol) this function works its way backwards to
        build the chain of the most important reactions (highest reaction rate for
        step) leading to said molecule in a given layer. It stops once it cannot find
        new molecules.'''
    k, m, r = max_re_one_step_layer(dat = dat, mol = mol, n = n, prod = True) # main production reaction for given species
    mol_checked = [mol]
    reaction_list = [mol + ':    ' +  r] # for easier readability and check, could have only r
    rate_list = [k]
    explored = False # this will turn True once the walker gets in a loop by checking the same molecules again
    while explored == False: 
        mol_next_layer = [] # let's work in layers, so need a list of molecules to check in each layer as going further back in the chain
        already_checked = 0
        for molec in m:
            if molec not in mol_checked and molec != 'M': # accounting for already checked molecules and third body
                k, mo, r = max_re_one_step_layer(dat = dat, mol = molec, n = n, prod = True)
                mol_checked.append(molec)
                reaction_list.append(molec + ':    ' + r) # for easier readability and check, could have only r
                rate_list.append(k)
                mol_next_layer += mo # these are lists so this should concatenate them
            elif molec in mol_checked:
                already_checked += 1
            
        if already_checked == len(m): # if every molecule to be checked has already been checked, exit by condition
            explored = True

        m = np.copy(mol_next_layer) # update m so next loop will check the next layer indeed, copy so python won't do weird stuff

    return reaction_list, rate_list

def all_layer_walker(dat, mol):
    ''' This function builds on the walker function and iterates through every layer.
        It returns a dictionary for each reaction found as keys, and their corresponding
        rate for each layer (in a list/array form) as values.
        If a reaction has not been listed in lower or higher layers the corresponding
        rates are considered NaN.'''
    react_dict = defaultdict(list) # new keys will automatically have empty lists
    for n in range(len(dat['atm']['pco'])):
        k, m, r = max_re_one_step_layer(dat = dat, mol = mol, n = n, prod = True) # main production reaction for given species
        react_dict[r].append(k)
        mol_checked = [mol]
        explored = False # this will turn True once the walker gets in a loop by checking the same molecules again
        rm = str(r) # basically copy, so it will not change
        while explored == False:
            mol_next_layer = [] # let's work in layers, so need a list of molecules to check in each layer as going further back in the chain
            already_checked = 0
            for molec in m:
                if molec not in mol_checked and molec != 'M': # accounting for already checked molecules and third body
                    k, mo, r = max_re_one_step_layer(dat = dat, mol = molec, n = n, prod = True)
                    mol_checked.append(molec)
                    if r not in react_dict.keys() and n >= 2: # checking for new reaction, in the first layer everything is new, so only need to check from 2nd layer
                        react_dict[r] = [np.nan for _ in range(len(react_dict[rm])-1)] # fill up missing values with nan
                    react_dict[r].append(k)
                    mol_next_layer += mo # these are lists so this should concatenate them
                elif molec in mol_checked:
                    already_checked += 1
                
            if already_checked == len(m): # if every molecule to be checked has already been checked, exit by condition
                explored = True

            m = np.copy(mol_next_layer) # update m so next loop will check the next layer indeed, copy so python won't do weird stuff
        
        for key in react_dict.keys(): # checking whether a reaction became not significant and if so, extend with nan
            if len(react_dict[key]) != n + 1:
                react_dict[key].append(np.nan)

    return react_dict
#%%
def plot_reactions(dat, mol, height = False, figname = None):
    ''' This function will execute all_layer_walker and then plot the reaction rates 
        for each reaction with pressure or height as y axis.
        
        need to make a colourscheme, see og plot scripts in vulcan'''
    y = []
    ylab = ''
    if height:
        y = dat['atm']['zco'][1:]/1e5 # height in km
        ylab = 'h [km]'
        yscale = 'lin'
    else:
        y = dat['atm']['pco']/1e6 # pressure in bar
        ylab = 'p [bar]'
        yscale = 'log'
    ps.reset_plt(23, 18)
    fig, ax = plt.subplots(tight_layout = True, figsize = (20,12))
    reaction_rates = all_layer_walker(dat, mol)
    colours = plt.cm.gist_ncar(np.linspace(0,1,len(reaction_rates.keys())))
    ls = ['-', '--', ':', '-.']
    for i,rea in enumerate(reaction_rates.keys()):
        if len(reaction_rates[rea]) == len(y): # until one small bug is found in all_layer_walker...
            ax.plot(reaction_rates[rea], y, label = rea, marker = 's', markersize = 1, color = colours[i], linestyle = ls[i%4])
    ax.set_xlabel(r'k [cm$^{-3}$ s$^{-1}$]')
    ax.set_ylabel(ylab)
    ax.set_xscale('log')
    ax.set_yscale(yscale)
    ax.legend(bbox_to_anchor = (1, 1), ncol = 2)
    ax.invert_yaxis()
    ax.set_xlim((1e-10,1e26))
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')

#%%
out_folder = '/scratch/s2555875/output/'
sim = out_folder + 'B_Pearce.vul' # hard coded for now, generalise later!!!-------------------------------

with open(sim, 'rb') as handle:
  data = pickle.load(handle)
#%%
R, K = walker(data, 'HCN', 0)

# %%
for i in range(len(R)):
    print('{}\t{:.1e}'.format(R[i],K[i]))
# %%
