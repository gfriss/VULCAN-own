#%%
import numpy as np
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
from petitRADTRANS.physics import guillot_global

atmosphere = Radtrans(line_species = ['H2O_HITEMP',
                                      'CH4',
                                      'CO2',
                                      'Na_allard',
                                      'K_allard',
                                      'N2'],
                      rayleigh_species = ['H2', 'He', 'N2'],
                      continuum_opacities = ['H2-H2', 'H2-He'],
                      wlen_bords_micron = [0.3, 15])

pressures = np.logspace(-1.301, 6, 120)
atmosphere.setup_opa_structure(pressures)

R_pl = nc.r_earth # maybe no need for these three
gravity = 980.
P0 = 1.

kappa_IR = 0.01
gamma = 0.4
T_int = 43.3
T_equ = 2534.5

temperature = guillot_global(pressures, kappa_IR, gamma, gravity, T_int, T_equ)

#%%
import matplotlib.pyplot as plt

plt.plot(temperature, pressures)
plt.yscale('log')
plt.show()