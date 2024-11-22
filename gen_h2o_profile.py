#%%
import numpy as np
import matplotlib.pyplot as plt
#import vulcan_cfg

pressure =  np.logspace(6,-2,120)
#pressure = pressure_dyne / 1.e6
#pressure = np.genfromtxt('atm_Earth_Jan_Kzz.txt', comments = '#', skip_header = 2)[:,0] / 1.e6
#pressure = np.genfromtxt('atm/atm_Earth_Jan_Kzz.txt', dtype = None, comments = '#', skip_header = 1, names = True)['Pressure']
# %%
def lin (x,a,b):
    return a*x + b

p_surface = pressure[0]
p_tropopause = 0.14*1e6 # in dyne
X_surface = 0.01
X_tropopause = 1.e-5
a_line = (X_tropopause - X_surface) / (p_tropopause - p_surface)
b_line = X_surface - a_line * p_surface

X_H2O = lin(pressure, a_line, b_line)
X_H2O[X_H2O < 1.e-6] = 1.e-6
#%%
# visual testing to match Pearce et al. (2022)
plt.plot(X_H2O, pressure/1e6)
plt.xlim(1.e-5,1.e-2)
plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()
# %%

new_file = '/home/s2555875/VULCAN-2/atm/mixing_table_archean.txt'
X_CO2 = 0.1
X_N2 = 0.1
X_CH4 = 5.e-3
X_O2 = 1.e-7
with open(new_file, 'w') as f:
    f.write('# (dyne/cm2)\nPressure  N2  CO2  CH4  O2  H2O\n')
    for i in range(len(pressure)):
        X_N2 = 1 - X_CO2 - X_CH4 - X_H2O[i] - X_O2
        f.write('{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\n'.format(pressure[i],X_N2,X_CO2,X_CH4,X_O2,X_H2O[i]))
        # simple test to make sure...
        if X_N2+X_CH4+X_CO2+X_H2O[i] != 1.:
            print('not perfect at level ' + str(i) + ' with sum of ' + str(X_N2+X_CH4+X_CO2+X_H2O[i]))
# %%
#table = np.genfromtxt(new_file, names=True, dtype=None, skip_header=1)
# %%
# generalised part with two lists as input (species, mixing ratio)
import sys

species = sys.argv[1].split(',')
mixing = sys.argv[2].split(',') # this may contain "fun" or something like this if it has to be calculated, e.g. not constant
out_file = sys.argv[3]
if 'func' in mixing:
    params = eval(sys.argv[4]) # need input format as '{"a":number,"b":number}'

def is_float(string):
    try: 
        float(string)
        return True
    except ValueError:
        return False

def gen_mix(x, p): # need to do general or something, this is only for linear
    a = (p['X_tropopause'] - p['X_surface']) / (p['p_tropopause'] - pressure[0])
    b = p['X_surface'] - a*pressure[0]
    X = a*x + b
    X[X < 0] = 0.0
    return X

nz = len(pressure)

ratios = {}
ratios['Pressure'] = pressure_dyne
for i,sp in enumerate(species):
    if mixing[i] == 'func':
        ratios[sp] = gen_mix(pressure, params) # this has to be changed according how the function works
    elif is_float(mixing[i]):
        ratios[sp] = np.ones(nz) * float(mixing[i])
    elif mixing[i] == 'last':
        rest = np.zeros(nz)
        for s in species:
            if s != sp:
                rest += ratios[s]
        ratios[sp] = np.ones(nz) - rest


import pandas as pd
df_ratios = pd.DataFrame.from_dict(ratios)

with open(out_file, 'w') as f:
    f.write('# (dyne/cm2)\n')
    f.write(df_ratios.to_string(header = True, index = False))
# %%
# calculating the total water mass used with a liner profile similar to that of Pearce et al. (2022)
# then distributing said water content with a constant mixing ratio in the atmosphere
from phy_const import kb, Navo
import numpy as np
from scipy.interpolate import interp1d
#%%
pT = np.genfromtxt('atm/TP_helios.txt', comments='#', skip_header=1, names=True, dtype = None)
pressure = pT['Pressure']
T = pT['Temp']
T_interp = interp1d(pressure, T)

pT_linwater = np.genfromtxt('atm/mixing_table_archean.txt', comments='#', skip_header=1, names=True, dtype = None)
pressure_linwater = pT_linwater['Pressure']
X_h2o_linwater = pT_linwater['H2O']
X_h2o_linwater = X_h2o_linwater[pressure_linwater >= np.min(pressure)]
pressure_linwater = pressure_linwater[pressure_linwater >= np.min(pressure)]
T_linwater = T_interp(pressure_linwater)
n0 = pressure_linwater / (kb*T_linwater) # gives total number density for each layer
total_n0 = np.sum(n0)
total_h2o = np.sum(n0*X_h2o_linwater)
X_h2o_const = total_h2o / total_n0
test_total_h2o = np.sum(n0*X_h2o_const)
test_X_h2o_const = test_total_h2o / total_n0
print('Calculated constant mixing ratio of water: {:.3e}'.format(X_h2o_const))
print('Testing conservation of water content:\n(original - constant mixing)/original = {} cm-3 '.format((total_h2o - test_total_h2o)/total_h2o))

# calculating constant mixing ratio of water for Archean following the
# total water column content of 20.2 kg m-2 by Wolf et al (2014)
N_h2o = 20.2 # kg/m2
N_h2o *= (1e-3 / 1e-4) # g/cm2
N_h2o *=  (Navo / 18) # 1/cm2, this is the integral of the waterprofile vertically so total_h2o previously
X_h2o_const_from_N = N_h2o/total_n0
print('Constant mixing ratio of water from column density: {:.3e}'.format(X_h2o_const_from_N))
# %%
X_co2 = 0.1
X_ch4 = 3e-3
X_o2 = 1e-7
X_n2 = 1 - X_co2 - X_ch4 - X_o2 - X_h2o_const
print('Constant mixing ratios using the integrated and divided water profile:')
print('N2: {:.3e}\nCO2: {:.3e}\nCH4: {:.3e}\nO2: {:.3e}\nH2O: {:.3e}'.format(X_n2, X_co2, X_ch4, X_o2, X_h2o_const))