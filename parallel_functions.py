#%%
import numpy as np
from astropy.constants import au,N_A,sigma_sb,L_sun,R_sun
from astropy import units as u
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
#%%
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
def gen_mixing(co2_mix, output):
    ''' Generates new initial mixing ratios and write them into a file compatable with VULCAN.
        It takes a base mixing ratio than changes the value for CO2, swapping it to CO
        which will change the C/O ratio.'''
    og_mixing = np.genfromtxt('atm/mixing_table_archean.txt', dtype = None, comments = '#', skip_header = 1, names = True)
    N2 = og_mixing['N2']
    H2O = og_mixing['H2O']
    CH4 = og_mixing['CH4']
    O2 = og_mixing['O2']
    CO2 = np.ones_like(N2) * co2_mix
    CO = np.ones_like(N2) * (np.max(co2_mix) - co2_mix)
    with open(output, 'w') as f:
        f.write('# (dyne/cm2)\nPressure  N2  CO2  CH4  O2  H2O  CO\n')
        for i in range(len(N2)):
            f.write('{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\n'.format(og_mixing['Pressure'][i],N2[i],CO2[i],CH4[i],O2[i],H2O[i],CO[i]))

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
        would not be too close to the edge (if needed).'''
    t_eff = df.loc[df.Name == name].T_eff.iloc[0]
    s_moist_gh = hz_inner_s_eff(t_eff)
    s_max_gh = hz_outer_s_eff(t_eff)
    s_eff = np.linspace(s_moist_gh/factor, s_max_gh*factor, n)
    return s_eff

# %%
def get_rad_prof(star):
    ''' To be used in run_parallel.py. Gives back the location of the needed stellar radiation profile
        file for a given star.'''
    rad_file = ''
    if star == 'EARLY_SUN':
        rad_file = 'atm/stellar_flux/Pearce_B_solar.txt'
    elif star == 'SUN':
        rad_file = 'atm/stellar_flux/Gueymard_solar.txt'
    else:
        rad_file = '/scratch/s2555875/stellar_flux/' + star.lower() + '.txt' # just to make sure it is lower case
    return rad_file
#%%
# function for local meteorite effect runs
def gen_mixing_local(h2_bar, output):
    ''' Generates new initial mixing ratios and write them into a file compatable with VULCAN.
        It takes a base mixing ratio than changes the value for CO2, swapping it to CO
        which will change the C/O ratio.
        
        .......'''
    og_mixing = np.genfromtxt('atm/mixing_table_archean.txt', dtype = None, comments = '#', skip_header = 1, names = True)
    N2 = og_mixing['N2'] # og was assumed that the partial pressure would be the same as mixing ratio
    H2O = og_mixing['H2O']
    CH4 = og_mixing['CH4']
    O2 = og_mixing['O2']
    CO2 = og_mixing['CO2']
    new_total = 1 + h2_bar # total pressure with new H2 amount
    N2 /= new_total
    H2O /= new_total
    CH4 /= new_total
    O2 /= new_total
    CO2 /= new_total
    H2 = np.ones_like(N2) * h2_bar / new_total
    with open(output, 'w') as f:
        f.write('# (dyne/cm2)\nPressure  N2  CO2  CH4  O2  H2O  H2\n')
        for i in range(len(N2)):
            f.write('{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\n'.format(og_mixing['Pressure'][i],N2[i],CO2[i],CH4[i],O2[i],H2O[i],H2[i]))