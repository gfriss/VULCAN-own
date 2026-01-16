#%%
import numpy as np
from astropy.constants import au,N_A,sigma_sb,L_sun,R_sun
from astropy import units as u
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import pickle
import os

scratch = '/scratch/s2555875' # place to store outputs

#%%
def get_str_number(n):
    ''' This function is used to get the number of digits in a number.'''
    if n < 10:
        return '0' + str(n)
    else:
        return str(n)
    
# C/O ratio using the 1 parameter study results
def get_element_number(species, element):
    ''' Gets the number of an element in a given species. Used to calculate C/O ratio later.'''
    char_list = list(species)
    ele_number = 0
    for i,char in enumerate(char_list):
        if i < len(char_list)-1 and char == element:
            if char_list[i+1].isnumeric():
                ele_number += 1*int(char_list[i+1])
            else:
                ele_number += 1
        if i == len(char_list)-1 and char == element:
            ele_number += 1
    return ele_number

def calc_C_to_O(dat, mixing_file):
    ''' It is generalised to get the species from the mixing ratio file used for the simulation.'''
    mix_data = np.genfromtxt(mixing_file, dtype = None, skip_header = 1, names = True)
    C_profile = np.zeros_like(dat['atm']['n_0'])
    O_profile = np.zeros_like(dat['atm']['n_0'])
    for name in mix_data.dtype.names:
        if name != 'Pressure' and 'C' in name:
            mul = get_element_number(name, 'C')
            C_profile += mul * dat['atm']['n_0'] * mix_data[name]
        if name != 'Pressure' and 'O' in name:
            mul = get_element_number(name, 'O')
            O_profile += mul * dat['atm']['n_0'] * mix_data[name]
    return trapezoid(C_profile, dat['atm']['zmco']) / trapezoid(O_profile, dat['atm']['zmco'])

def get_C_to_O_conjoint(nsim, network, nowash, version):
    ''' Using the above function, calculate the C/O ratio for all runs in the conjoint study.'''
    C_to_O = []
    for j in range(nsim['CtoO']):
        sim_CtoO = 'sim_{:02d}_CtoO'.format(j)
        mixing_file = os.path.join(scratch, 'mixing_files', sim_CtoO + version + 'mixing.txt')
        with open(os.path.join(scratch, 'output', sim_CtoO + network + nowash + version + '.vul'), 'rb') as handle:
            data_CtoO = pickle.load(handle)
        C_to_O.append(round(calc_C_to_O(data_CtoO, mixing_file), 4)) # round to 4 decimal places to avoid issues with different machines
    return C_to_O

# BC from meteorite bombardment rate:
def bar_to_number(x): 
    ''' Converts partial pressure of a species to number density.'''
    return x * 2300 * N_A*u.mol / 100 # 100 bar is 2300 mol/cm-2 according to Zahnler et al. (2020)

def gpergyr_to_gpers(x):
    ''' g/yr to g/s.'''
    return x / 3.154e16

def bombardment(rate, mass, prod):
    ''' Function to calculate the surface emission due to meteoritic bombardment using Zahnler et al. (2020)
        equilibrium production (Table1 last rows). Rate is in g/Gyr, mass is in g and prod is a dict
        containing the molecules and their eq. production in bars.'''
    bc = {}
    for sp in prod.keys():
        nd_per_mass = bar_to_number(prod[sp]) / mass # 1/(cm2 g)
        surface_flux = nd_per_mass * gpergyr_to_gpers(rate)
        bc[sp] = str(surface_flux)
    return bc

def dict_to_input(d):
    ''' Converts a dictionary into input for a python script, namely gen_cfg.py that take
        the input as key1,key2,etc. val1,val2,etc.'''
    k = ','.join(d.keys())
    v = []
    for val in d.values():
        v.append(val)
    val_str = ','.join(v)
    return k, val_str

# initial mixing ratio for C/O tests
def gen_mixing(new_mix, output, species = 'CO2', version = ''):
    ''' Generates new initial mixing ratios and write them into a file compatable with VULCAN.
        It takes a base mixing ratio than changes the value for CO2, swapping it to CO
        which will change the C/O ratio, or changing CH4 (potentially adding something else to balance?).'''
    og_mixing = np.genfromtxt('atm/mixing_table_archean{}.txt'.format(version), dtype = None, comments = '#', skip_header = 1, names = True)
    N2 = og_mixing['N2']
    H2O = og_mixing['H2O']
    O2 = og_mixing['O2']
    H2 = np.zeros_like(N2)
    CO = np.zeros_like(N2)
    if 'CO' in og_mixing.dtype.names:
        CO = og_mixing['CO']
    if species == 'CO2': # change Co2 and add CO, keep CH4 the same
        CO2 = np.ones_like(N2) * new_mix
        CO += np.ones_like(N2) * (0.1 - new_mix)
        CH4 = og_mixing['CH4']
    elif species == 'CH4': # change CH4 and keep CO2 the same, set CO to 0 for universal use of writing part later
        CH4 = np.ones_like(N2) * new_mix
        CO2 = og_mixing['CO2']
    elif species == 'CH4-H2': # change both CH4 and H2 to somewhat keep H2 and redox balance, set CO to 0
        CH4 = np.ones_like(N2) * new_mix
        CO2 = og_mixing['CO2']
        H2 = np.ones_like(N2) * (5.e-3 - new_mix) # keeping total of 5000 ppmv for CH4 + H2
    elif species == 'CH4-balanced': # change CH4 but introduce CO and decrease CO2 (to keep C/O ratio) and H2O (to keep H content, somewhat)
        with open(os.path.join(scratch, 'output/archean_ncho.vul'), 'rb') as handle:
            d = pickle.load(handle)
        original_ctoo = calc_C_to_O(d, 'atm/mixing_table_archean.txt')
        n0 = d['atm']['n_0']        
        CH4 = np.ones_like(N2) * new_mix
        CO = np.ones_like(N2) * (original_ctoo*(np.sum(n0*2*(og_mixing['CO2']+og_mixing['O2'])) + np.sum(n0*og_mixing['H2O'])) - np.sum(n0*(og_mixing['CO2']+CH4)) ) / (2*original_ctoo*np.sum(n0))
        CO2 = og_mixing['CO2'] - CO
        N2 = np.ones_like(N2) - (CO2 + CH4 + O2 + H2O + CO + H2)
        #H2 = np.ones_like(N2) - (N2 + CO2 + CH4 + O2 + H2O)
    elif species == 'CH4-updated':
        CH4 = np.ones_like(N2) * new_mix
        CO = 100. * CH4
        CO2 = og_mixing['CO2']
        N2 = np.ones_like(N2) - (CO2 + CH4 + O2 + H2O + CO + H2)
    with open(output, 'w') as f:
        f.write('# (dyne/cm2)\nPressure  N2  CO2  CH4  O2  H2O  CO  H2\n')
        for i in range(len(N2)):
            f.write('{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\n'.format(og_mixing['Pressure'][i],N2[i],CO2[i],CH4[i],O2[i],H2O[i],CO[i], H2[i]))
            
# rainout rate calculation
def rainout(dat, rain_spec = 'HCN_rain', g_per_mol = 27):
    ''' Calculates the rainout rate of the given species and returns it with units of kg/m2/yr.'''
    #rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = trapezoid(dat['variable']['y_rain'][rain_spec], dat['atm']['zmco']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = rain_rate * 5.237e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

# functions for different stellar type tests
def read_stellar_data(file):
    ''' Stellar data is saved in a csv file and has str and float values. This read in function
        makes sure that every read in value has the correct type.'''
    types = defaultdict(lambda: 'float', Name = 'str', Type = 'str', T_eff_source = 'str', L_R_source = 'str', Dist_source = 'str')
    return pd.read_csv(file, dtype = types)

def calc_T_eq(L_star, albedo = 0.06, a = au): 
    ''' Calculates standard sim's T_eq.'''
    return np.power((1-albedo)*L_star / (16*np.pi*sigma_sb*(a**2)), 1/4) #something is wrong gives too high value...proceed with R and T eq. for now

def semi_major_axis(T_star, R_star, T_eq = 1, albedo = 0.06):
    ''' Calculates the semi major axis given the stellar effective temperature and the planets equilibrium temperature.'''
    a = (np.power(T_star, 2)/np.power(T_eq, 2)) * np.sqrt(1-albedo) * (R_star/2)
    return a / au # m to AU
#%%
def test_semi_major_axis(df):
    ''' Test function to calculate the semi major axis and recieved flux for all planets in 
        simulations with the equilibrium temperature approach.'''
    T_sun_ee = np.power(L_sun / (4*np.pi*(R_sun**2)*sigma_sb), 0.25)
    T_eq_ee = T_sun_ee * np.power(1-0.06, 0.25) * np.sqrt(R_sun/(2*au))
    a_list = []
    s_eff_list = []
    for t,r,l in zip(df.T_eff,df.R,df.L_log):
        d = semi_major_axis(t*u.K, r*R_sun, T_eq=T_eq_ee)
        #print(t*u.K, '->', d*u.au)
        #print('S_eff: ', t*u.K, '->', np.power( np.power(10, l) / d, 2) )
        a_list.append(d)
        s_eff_list.append(np.power( np.power(10, l) / d, 2))
    
    return a_list, s_eff_list

def hz_inner_s_eff(temp_star_eff):
    ''' Calculates the inner edge of the habitable zone (moist greenhouse limit) described in Kopprapau et al. (2013).'''
    temp_star = temp_star_eff - 5780
    return 1.014 + 8.1774e-5*temp_star + 1.7063e-9*(temp_star**2) - 4.3241e-12*(temp_star**3) - 6.6462e-16*(temp_star**4)

def hz_outer_s_eff(temp_star_eff):
    '''  Calculates the outer edge of the habitable zone (mximum greenhouse limit) described in Kopprapau et al. (2013).'''
    temp_star = temp_star_eff - 5780
    return 0.3438 + 5.8942e-5*temp_star + 1.6558e-9*(temp_star**2) - 3.0045e-12*(temp_star**3) - 5.2983e-16*(temp_star**4)

def plot_hz_for_test(df, s_eff, figname = None):
    ''' Following Kopprapau et al. (2013) it ests whether the given approach puts the planets into similar part of 
        habitable zone around various stars using the effective flux (normalised by current Earth value).'''
    T_star_eff = np.linspace(2500,7000,42)
    s_eff_inner = hz_inner_s_eff(T_star_eff)
    s_eff_outer = hz_outer_s_eff(T_star_eff)
    fig, ax = plt.subplots(tight_layout = True)
    ax.plot(s_eff_inner, T_star_eff, 'r', label = 'Moist greenhouse')
    ax.plot(s_eff_outer, T_star_eff, 'b--', label = 'Maximum greenhouse')
    ax.plot(s_eff, df.T_eff, 'go', linestyle = '', label = 'Planets')
    ax.invert_xaxis()
    ax.set_xlabel(r'Effective flux incident on planet (S/S$_0$)')
    ax.set_ylabel(r'T$_{eff}$ [K]')
    ax.legend()
    if figname != None:
        fig.savefig('/scratch/s2555875/plot/' + figname)

def test_semi_major_from_S_eff(df, f = 0.744):
    ''' Test function to calculate the semi major axis and recieved flux for all planets in 
        simulations with the effective flux approach. This approach takes the ratio of the
        effective flux of the base reaction and the moist greenhouse limit (f) and places
        every planet such that they all share this ratio which keeps them in the same
        part of the habitable zone.'''
    a_list, s_eff_list = [], []
    for l,t in zip(df.L_log,df.T_eff):
        s_moist_gh = hz_inner_s_eff(t)
        s_eff = f * s_moist_gh
        s_eff_list.append(s_eff)
        d = np.sqrt( np.power(10, l) / s_eff )
        a_list.append(d)
    return a_list, s_eff_list

def semi_major_from_S_eff(df, name, f = 0.744):
    ''' Same as before, but only for one planet. To be used in run_parallel.py.'''
    t_eff = df.loc[df.Name == name].T_eff.iloc[0]
    Llog = df.loc[df.Name == name].L_log.iloc[0]
    s_moist_gh = hz_inner_s_eff(t_eff)
    s_eff = f * s_moist_gh
    a = np.sqrt( np.power(10, Llog) / s_eff )
    return a

def semi_major_list_from_Seff(df, name, n, factor = 1):
    ''' Gives back a range of distances in the habitable zone. Factor changes the vaules so
        would not be too close to the edge (if needed).'''
    t_eff = df.loc[df.Name == name].T_eff.iloc[0]
    Llog = df.loc[df.Name == name].L_log.iloc[0]
    s_moist_gh = hz_inner_s_eff(t_eff)
    s_max_gh = hz_outer_s_eff(t_eff)
    s_eff = np.linspace(s_moist_gh/factor, s_max_gh*factor, n)
    a = np.sqrt( np.power(10, Llog) / s_eff )
    return a

def Seff_list(df, name, n, factor = 1):
    ''' Gives back a range of effective radiations in the habitable zone. Factor changes the vaules so
        would not be too close to the edge (if needed, should be >=1).'''
    t_eff = df.loc[df.Name == name].T_eff.iloc[0]
    s_moist_gh = hz_inner_s_eff(t_eff)
    s_max_gh = hz_outer_s_eff(t_eff)
    s_eff = np.linspace(s_moist_gh*factor, s_max_gh/factor, n)
    return s_eff

# %%
def get_rad_prof(star):
    ''' Gives back the location of the needed stellar radiation profile file for a given star.'''
    rad_file = ''
    if star == 'SUN':
        rad_file = '/home/s2555875/VULCAN-2/atm/stellar_flux/Gueymard_solar.txt'
    else:
        rad_file = '/scratch/s2555875/stellar_flux/' + star.lower() + '.txt' # just to make sure it is lower case
    return rad_file
#%%
def get_sim_not_done(file, star_df, param_dict, C_to_O, nsim, nsim_total):
    sim_done = []
    if os.path.isfile(os.path.join(scratch, file)):
        with open(os.path.join(scratch, file), 'r') as f:
            for line in f:
                if line[0].isnumeric():
                    bits = line.split()[:3] # keeping Teff, dist and CtoO
                    star_name = star_df.loc[star_df.T_eff == float(bits[0]), 'Name'].iloc[0]
                    i_dist = param_dict[star_name].index(float(bits[1]))
                    i_CtoO = C_to_O.index(float(bits[2]))
                    sim_done.append((star_name, i_dist, i_CtoO))
    sim_not_done = []
    for i in range(nsim_total):
        i_star = i//(nsim['dist']*nsim['CtoO']) # which star (changes after looped through all distance and C/O possibilities)
        star_name = star_df.Name.loc[i_star]
        i_dist = (i//nsim['CtoO'])%nsim['dist'] # which distance (changes after looped through all C/O possibilities and then restarts when all distances are done)
        i_CtoO = i%nsim['CtoO'] # which C/O (simply loops through the C/O possibilities)
        if (star_name, i_dist, i_CtoO) not in sim_done:
            sim_not_done.append((star_name, i_dist, i_CtoO))
    return sorted(sim_not_done)