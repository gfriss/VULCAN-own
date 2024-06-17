#%%
import numpy as np
import pickle

def calc_rate(dat, spec_list, re_id, n):
    rate = dat['variable']['k'][re_id][n]
    for sp in spec_list:
        if sp != 'M':
            rate *= dat['variable']['y'][n, dat['variable']['species'].index(sp)]
    return rate

def get_species(eq_side):
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
    k, m, r = max_re_one_step_layer(dat = dat, mol = mol, n = n, prod = True) # main production reaction for given species
    mol_checked = [mol]
    reaction_list = [mol + ':    ' +  r]
    rate_list = [k]
    explored = False # this will turn True once the walker gets in a loop by checking the same molecules again
    while explored == False: 
        mol_next_layer = [] # let's work in layers, so need a list of molecules to check in each layer as going further back in the chain
        already_checked = 0
        for molec in m:
            if molec not in mol_checked and molec != 'M': # accounting for already checked molecules and third body
                k, mo, r = max_re_one_step_layer(dat = dat, mol = molec, n = n, prod = True)
                mol_checked.append(molec)
                reaction_list.append(molec + ':    ' + r)
                rate_list.append(k)
                mol_next_layer += mo # these are lists so this should concatenate them
            elif molec in mol_checked:
                already_checked += 1
            
        if already_checked == len(m): # if every molecule to be checked has already been checked, exit by condition
            explored = True

        m = np.copy(mol_next_layer) # update m so next loop will check the next layer indeed, copy so python won't do weird stuff

    return reaction_list, rate_list
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
