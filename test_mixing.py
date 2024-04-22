#%%
import numpy as np
import matplotlib.pyplot as plt

og_mixing = np.genfromtxt('atm/mixing_Pearce_B.txt', dtype = None, comments = '#', skip_header = 1, names = True)
#og_CtoO = np.sum(og_mixing['CO2'] + og_mixing['CH4']) / np.sum(2*og_mixing['CO2'] + og_mixing['H2O'])

N2 = og_mixing['N2']
H2O = og_mixing['H2O']
CH4 = og_mixing['CH4']
CO2_mix = np.linspace(0,0.9, n, endpoint = True)
#%%
c_to_o = []
for mix in CO2_mix:
    CO2 = np.ones_like(N2) * mix
    CO = np.ones_like(N2) * (0.9-mix) # max 90%
    c_to_o.append(np.sum(CO2 + CH4 + CO) / np.sum(CO + H2O + 2*CO2))

plt.plot(CO2_mix, c_to_o)
# %%
