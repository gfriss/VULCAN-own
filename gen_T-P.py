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

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#%%
beta_1 = 0.5
beta_2 = 0.5

def T_from_p(p, p0, p1, p2, p3, T0, T1, T2, T3):
    alpha_2 = np.log(p3/p2) / (np.power(T3-T2, beta_2))
    alpha_1 = (np.log(p2/p0) + alpha_2*np.power(T1-T2, beta_2)) / (np.power(T1-T0, beta_1))
    print(alpha_1)
    print(alpha_2)
    
    T = np.zeros_like(p)
    T[p < p1] = T0 + np.power((1/alpha_1)*np.log(p[p < p1]/p0), 1/beta_1)
    T[p >= p1] = T2 + np.power((1/alpha_2)*np.log(p[p >= p1]/p2), 1/beta_2)
    return T

pressures = np.logspace(-7.301, 0, 71) # converted to bars already
#%%
temperature = T_from_p(pressures, pressures[0], 6e-1, 6.5e-1, pressures[-1], 214., 215.1, 215., 323.) # for ox it was 6e-2, 6e-2 and 351

plt.plot(temperature, pressures)
plt.yscale('log')
plt.gca().invert_yaxis()    
# %%
vulcan_atm = np.genfromtxt('atm/atm_Earth_Jan_Kzz.txt', dtype = None, comments = '#', skip_header = 1, names = True)
Kzz = vulcan_atm['Kzz']
n = len(Kzz)

pressures = vulcan_atm['Pressure'] / 1e6
temperature = 109 * np.power(pressures, 4) + 214.

plt.plot(temperature, pressures)
plt.yscale('log')
plt.gca().invert_yaxis()
#%%
sclae_Kzz = 6.3
with open('../VULCAN-master/atm/T-P_adiabat.txt', 'w') as f:
    f.write('# (dyne/cm2) (K)     (cm2/s)\nPressure\tTemp\tKzz\n')
    for i in range(len(temperature)):
        new_Kzz = Kzz[i]*sclae_Kzz
        if i < 10 and new_Kzz > 6.3e5: # 10 is arbitrary, but essentially making sure that first few behaves as in Pearce et al.
            new_Kzz = 6.3e5
        f.write('{:.3e}\t{:.1f}\t{:.3e}\n'.format(pressures[i]*1e6, temperature[i], new_Kzz))

# %%
# fitting the parametric T-P model
p_0 = pressures[0]
p_3 = pressures[-1] 
       
def Tp_fit(x, p1, p2, T0, T2, a1, a2):
    T = np.zeros_like(x)
    T[x < p1] = T0 + np.power((1/a1)*np.log(x[x < p1]/p_0), 1/beta_1)
    T[x >= p1] = T2 + np.power((1/a2)*np.log(x[x >= p1]/p2), 1/beta_2)
    return T
#%%
pres = pressures[::3]
temp = temperature[::3]
rng = np.random.default_rng()
temp += 0.2 * rng.normal(size=pres.size)

popt, pcov = curve_fit(Tp_fit, pres, temp, p0 = [5e-1, 8e-1,210.,215., 14, 0.01])
perr = np.diag(pcov)

print('p1, p2, T0, T2, a1, a2\n')
print(popt)
print('\n+-\n')
print(perr)

plt.plot(Tp_fit(pressures, *popt), pressures)
plt.yscale('log')
plt.gca().invert_yaxis()  
# %%
# adiabatic version
pressures = np.logspace(-7.301, 0, 71)
temperature = 137 * np.power(pressures, 1.2) + 214.

plt.plot(temperature, pressures)
plt.axhline(y=pressures[-11], linestyle='--')
plt.axhline(y=pressures[-12], linestyle='-.')
plt.axhline(y=pressures[-13])
plt.yscale('log')
plt.gca().invert_yaxis()

# %%
# testing the ddy diffucion profile
import numpy as np
import matplotlib.pyplot as plt

vulcan_atm = np.genfromtxt('atm/atm_Earth_Jan_Kzz.txt', dtype = None, comments = '#', skip_header = 1, names = True)
Kzz = vulcan_atm['Kzz']
n = len(Kzz)

pressures = vulcan_atm['Pressure'] / 1e6
press_new = np.logspace(-7.301, 0, n)

plt.plot(Kzz, pressures, label = 'VULCAN')
plt.plot(Kzz, press_new[::-1], '--', label = 'new')
plt.yscale('log')
plt.gca().invert_yaxis()
# %%
# scaling the current Earth eddy diffusion profile according to Hu et al. (2012)
# and using the T-P profile from Pearce et al. (2022)
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
scale = 6.3 # for H2 dominated atmosphere, 0.68 for CO2 dominated

vulcan_atm = np.genfromtxt('atm/atm_Earth_Jan_Kzz.txt', dtype = None, comments = '#', skip_header = 1, names = True)
Kzz_vulcan = vulcan_atm['Kzz']
pressure_vulcan = vulcan_atm['Pressure']
pearce_atm = np.genfromtxt('atm/T-P-Kzz_Pearce_A.txt', dtype = None, comments = '#', skip_header = 1, names = True) # already has vulcan format, but eddy diffusion is weird
pressure_pearce = pearce_atm['Pressure']
temperature_pearce = pearce_atm['Temp']

kzz_fun = interp1d(pressure_vulcan, Kzz_vulcan*scale)
Kzz_scaled_overlap = kzz_fun(pressure_pearce[(pressure_pearce > np.min(pressure_vulcan))*(pressure_pearce < np.max(pressure_vulcan))])
Kzz_scaled_high_pressure = np.ones_like(pressure_pearce[pressure_pearce > np.max(pressure_vulcan)]) * Kzz_scaled_overlap[0] # by structure it starts at ground so max pressure
Kzz_scaled_low_pressure = np.ones_like(pressure_pearce[pressure_pearce < np.min(pressure_vulcan)]) * Kzz_scaled_overlap[-1] # by structure it ends at TOA so min pressure
Kzz_scaled = np.concatenate((np.concatenate((Kzz_scaled_high_pressure, Kzz_scaled_overlap)), Kzz_scaled_low_pressure))

# visual test
plt.plot(Kzz_scaled, pressure_pearce)
plt.gca().invert_yaxis()
plt.xscale('log')
plt.yscale('log')

with open('atm/T-P-Kzz_scaled_A.txt', 'w') as f:
    f.write('# (dyne/cm2) (K)     (cm2/s)\nPressure\tTemp\tKzz\n')
    for i in range(len(temperature_pearce)):
        f.write('{:.3e}\t{:.1f}\t{:.3e}\n'.format(pressure_pearce[i], temperature_pearce[i], Kzz_scaled[i]))
# %%
