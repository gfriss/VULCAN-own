#%%
import numpy as np
import matplotlib.pyplot as plt
import vulcan_cfg

pressure_dyne =  np.logspace(np.log10(vulcan_cfg.P_b),np.log10(vulcan_cfg.P_t),vulcan_cfg.nz)
pressure = pressure_dyne / 1.e6
#pressure = np.genfromtxt('atm_Earth_Jan_Kzz.txt', comments = '#', skip_header = 2)[:,0] / 1.e6

# %%
def lin (x,a,b):
    return a*x + b

p_surface = pressure[0]
p_tropopause = 0.14
X_surface = 0.01
X_tropopause = 1.e-5
a_line = (X_tropopause - X_surface) / (p_tropopause - p_surface)
b_line = X_surface - a_line * p_surface

X_H2O = lin(pressure, a_line, b_line)
X_H2O[X_H2O < 0.] = 0.

#%%
# visual testing to match Pearce et al. (2022)
#plt.plot(X_H2O, pressure)
#plt.xlim(1.e-5,1.e-2)
#plt.yscale('log')
#plt.xscale('log')
#plt.gca().invert_yaxis()
# %%

new_file = '/home/s2555875/VULCAN-2/atm/mixing_table.txt'

X_CO2 = 0.6
X_N2 = 0.1
X_CH4 = 10.e-6
X_Ar = 9.34e-3
X_SO2 = 2e-10
with open(new_file, 'w') as f:
    f.write('# (dyne/cm2)\nPressure  CO2  N2  CH4  H2O  Ar  SO2\n')
    for i in range(len(pressure)):
        #X_N2 = (1 - X_CH4 - X_H2O[i] - X_Ar - X_SO2) / 10 # 1x+9x for N2 and H2
        #X_CO2 = 9 * X_N2
        X_N2 = 1 - X_CO2 - X_CH4 - X_Ar - X_SO2 - X_H2O[i]
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pressure_dyne[i],X_CO2,X_N2,X_CH4,X_H2O[i],X_Ar,X_SO2))
# %%
#table = np.genfromtxt(new_file, names=True, dtype=None, skip_header=1)

#n_H = 4 * np.sum(table['CH4']) + 2 * np.sum(table['H2O'])
#n_O = 2 * np.sum(table['CO2']) + 2 * np.sum(table['SO2'])
#n_C = np.sum(table['CO2']) + np.sum(table['CH4'])
#n_N = 2 * np.sum(table['N2'])
#n_S = np.sum(table['SO2'])

#f_O = n_O / n_H
#f_C = n_C / n_H
#f_N = n_N / n_H
#f_S = n_S / n_H
# %%
#print('f_O = {:.4e}'.format(f_O))
#print('f_C = {:.4e}'.format(f_C))
#print('f_N = {:.4e}'.format(f_N))
#print('f_S = {:.4e}'.format(f_S))
# %%
