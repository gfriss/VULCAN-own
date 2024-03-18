#%%
import numpy as np
import matplotlib.pyplot as plt
import vulcan_cfg

pressure =  np.logspace(np.log10(vulcan_cfg.P_b),np.log10(vulcan_cfg.P_t),vulcan_cfg.nz)
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
X_H2O[X_H2O < 0.] = 0.
#%%
# visual testing to match Pearce et al. (2022)
plt.plot(X_H2O, pressure/1e6)
plt.xlim(1.e-5,1.e-2)
plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()
# %%

new_file = '/home/s2555875/VULCAN-master/atm/mixing_table.txt'
X_CO2 = 0.9
X_N2 = 0.1
X_CH4 = 2.e-6
X_Ar = 9.34e-3
X_SO2 = 2e-10
with open(new_file, 'w') as f:
    f.write('# (dyne/cm2)\nPressure  CO2  N2  CH4  H2O\n')#  SO2\n')#  Ar\n')
    for i in range(len(pressure)):
        #X_N2 = (1 - X_CH4 - X_H2O[i] - X_Ar - X_SO2) / 10 # 1x+9x for N2 and H2
        #X_H2 = 9 * X_N2
        X_N2 = 1 - X_CO2 - X_CH4 - X_H2O[i]# - X_SO2# - X_Ar
        f.write('{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\t{:.3E}\n'.format(pressure[i],X_CO2,X_N2,X_CH4,X_H2O[i]))#)),X_SO2))#,X_Ar))
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
