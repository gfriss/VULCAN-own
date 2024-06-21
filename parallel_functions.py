import numpy as np
from astropy.constants import au,N_A

# BC from meteorite bombardment rate:
def bar_to_number(x): # 100 bar is 2300 mol/cm-2 according to Zahnler et al. (2020)
    return x * 2300 * N_A / 100

def gpergyr_to_gpers(x):
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
    k = ','.join(d.keys())
    v = []
    for val in d.values():
        v.append(val)
    val_str = ','.join(v)
    return k, val_str

def gen_mixing(co2_mix, output):
    og_mixing = np.genfromtxt('atm/mixing_Pearce_B.txt', dtype = None, comments = '#', skip_header = 1, names = True)
    N2 = og_mixing['N2']
    H2O = og_mixing['H2O']
    CH4 = og_mixing['CH4']
    CO2 = np.ones_like(N2) * co2_mix
    CO = np.ones_like(N2) * (0.9 - co2_mix) # max 90 %
    with open(output, 'w') as f:
        f.write('# (dyne/cm2)\nPressure  CO2  CO  N2  CH4  H2O\n')
        for i in range(len(N2)):
            f.write('{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\n'.format(og_mixing['Pressure'][i],CO2[i],CO[i],N2[i],CH4[i],H2O[i]))

def semi_major_axis(T_star, R_star, T_eq = 1, albedo = 0.06):
    a = (np.power(T_star, 2)/np.power(T_eq, 2)) * np.sqrt(1-albedo) * (R_star/2) # cm
    return a / (au/1000.) # cm to AU