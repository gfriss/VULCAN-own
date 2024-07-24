#%%
import numpy as np
from petitRADTRANS.radtrans import Radtrans
#from petitRADTRANS import nat_cst as nc
#from petitRADTRANS.physics import guillot_global
import petitRADTRANS.physical_constants as cst
from petitRADTRANS.spectral_model import SpectralModel
import matplotlib.pyplot as plt
import petitRADTRANS.chemistry.utils as cu
from scipy.interpolate import interp1d
import pickle
#%%
mix_name = '/scratch/s2555875/HELIOS/input/vmr_mix_earth.txt'
mixing_data = np.genfromtxt(mix_name, dtype = None, comments = '#', skip_header = 1, names = True)
mixing = {}
pres = np.logspace(np.log10(5e-8), np.log10(1), 100, endpoint = True)
for name in mixing_data.dtype.names:
    if name != 'Pressure':
        int_mix = interp1d(mixing_data['Pressure']/1e6, mixing_data[name], bounds_error = False, fill_value = 'extrapolate')
        mixing[name] = int_mix(pres)

mass_ratios = cu.volume_mixing_ratios2mass_fractions(mixing)

earth_vulcan = '/scratch/s2555875/output/Earth.vul'
with open(earth_vulcan, 'rb') as handle:
    earth_data = pickle.load(handle)
hel_tp = np.genfromtxt('/scratch/s2555875/HELIOS/output/earth/earth_tp.dat', dtype = None, names = True, skip_header = 1,usecols=(1,2))
#%%
# current Earth comparison
spectral_model = SpectralModel(
    # Radtrans parameters
    pressures=pres,
    line_species=[
        'H2O',
        'CH4',
        'CO2',
        'O2',
        'O3'
    ],
    rayleigh_species=['N2', 'O2'],
    gas_continuum_contributors=['N2-N2', 'N2-O2', 'O2-O2', 'H2O-N2', 'H2O-H2O'],
    wavelength_boundaries=[0.05, 100],
    scattering_in_emission=True,  # replace do_scat_emis from pRT2
    # SpectralModel parameters
    # Planet parameters
    planet_radius=1 * cst.r_earth,
    reference_gravity=980,
    reference_pressure=1,
    # Star, system, orbit
    is_observed=False,  # return the flux observed at system_distance
    is_around_star=True,  # if True, calculate a PHOENIX stellar spectrum and add it to the emission spectrum
    #system_distance=10 * cst.s_cst.light_year * 1e2,  # m to cm, used to scale the spectrum
    star_effective_temperature=5780,  # used to get the PHOENIX stellar spectrum model
    star_radius=1 * cst.r_sun,  # used to get the star flux irradiating the planet
    orbit_semi_major_axis=1 * cst.au,  # used to get the star flux irradiating the planet
    # Temperature profile parameters
    temperature_profile_mode='guillot',
    temperature = 274.,
    intrinsic_temperature=43.3,
    guillot_temperature_profile_gamma=0.3,
    guillot_temperature_profile_kappa_ir_z0=0.00196,
    #metallicity = 1,
    #co_ratio = 0.47,
    # Mass fractions
    #use_equilibrium_chemistry=False,
    imposed_mass_fractions=mass_ratios,
    filling_species={  # automatically fill the atmosphere with H2 and He, such that the sum of MMRs is equal to 1 and H2/He = 37/12, H2/Ne = 37/0.06
        'H2': 10,
        'He': 2,
        'Ne': 0.1,
        'Ar': 0.1
    }
)

wavelengths, flux = spectral_model.calculate_spectrum(
    mode='emission',
    update_parameters=True  # this will build notably the temperature and mass fractions profile
)

fig, ax = plt.subplots(figsize = (10,6))
ax.plot(earth_data['atm']['Tco'], earth_data['atm']['pco']/1e6, label = 'VULCAN')
ax.plot(hel_tp['tempK'], hel_tp['press106bar']/1e6, label = 'HELIOS', linestyle = '--')
ax.plot(spectral_model.temperatures, spectral_model.pressures * 1e-6, label = 'petitRADTRANS', linestyle = ':')
ax.set_yscale('log')
ax.set_ylim([1, 5e-8])
ax.set_xlabel('T [K]')
ax.set_ylabel('P [bar]')
ax.legend()
fig.savefig('/scratch/s2555875/plot/current_earth_TP_comp.png')
#%%
# Early Earth comparison
tp_pearce = np.genfromtxt('atm/T-P-Kzz_Pearce_B.txt', dtype = None, names = True, skip_header = 1, usecols=(0,1))
hel_tp_early = np.genfromtxt('/scratch/s2555875/HELIOS/output/3/3_tp.dat', dtype = None, names = True, skip_header = 1,usecols=(1,2))
pres_early = np.logspace(-8, np.log10(2), 100, endpoint = True)
mix_name_early = '/scratch/s2555875/HELIOS/input/vmr_mix.txt'
mixing_data_early = np.genfromtxt(mix_name, dtype = None, comments = '#', skip_header = 1, names = True)
mixing_early = {}

for name in mixing_data_early.dtype.names:
    if name != 'Pressure':
        int_mix_early = interp1d(mixing_data_early['Pressure']/1e6, mixing_data_early[name], bounds_error = False, fill_value = 'extrapolate')
        mixing_early[name] = int_mix_early(pres)

mass_ratios_early = cu.volume_mixing_ratios2mass_fractions(mixing_early)
#%%
spectral_model = SpectralModel(
    # Radtrans parameters
    pressures=pres_early,
    line_species=[
        'H2O',
        'CH4',
        'CO2'
    ],
    rayleigh_species=['N2', 'CO2'],
    gas_continuum_contributors=['CO2-CO2', 'N2-N2', 'H2O-H2O', 'H2O-N2'],
    wavelength_boundaries=[0.05, 100],
    scattering_in_emission=True,  # replace do_scat_emis from pRT2
    # SpectralModel parameters
    # Planet parameters
    planet_radius=1 * cst.r_earth,
    reference_gravity=980,
    reference_pressure=2,
    # Star, system, orbit
    is_observed=False,  # return the flux observed at system_distance
    is_around_star=True,  # if True, calculate a PHOENIX stellar spectrum and add it to the emission spectrum
    #system_distance=10 * cst.s_cst.light_year * 1e2,  # m to cm, used to scale the spectrum
    star_effective_temperature=5332,  # used to get the PHOENIX stellar spectrum model
    star_radius=1 * cst.r_sun,  # used to get the star flux irradiating the planet
    orbit_semi_major_axis=1 * cst.au,  # used to get the star flux irradiating the planet
    # Temperature profile parameters
    temperature_profile_mode='guillot',
    temperature = 253.15,
    intrinsic_temperature=43.3,
    guillot_temperature_profile_gamma=0.0005,
    guillot_temperature_profile_kappa_ir_z0=0.0015,
    # Mass fractions
    imposed_mass_fractions=mass_ratios_early,
    filling_species={  # automatically fill the atmosphere with H2 and He, such that the sum of MMRs is equal to 1 and H2/He = 37/12, H2/Ne = 37/0.06
        'H2': 10,
        'He': 2,
        'Ne': 0.1,
        'Ar': 0.1
    }
)

wavelengths, flux = spectral_model.calculate_spectrum(
    mode='emission',
    update_parameters=True  # this will build notably the temperature and mass fractions profile
)

fig, ax = plt.subplots(figsize = (10,6))
ax.plot(tp_pearce['Temp'], tp_pearce['Pressure']/1e6, label = 'Pearce et al. (2022)')
ax.plot(hel_tp_early['tempK'], hel_tp_early['press106bar']/1e6, label = 'HELIOS', linestyle = '--')
ax.plot(spectral_model.temperatures, spectral_model.pressures * 1e-6, label = 'petitRADTRANS', linestyle = ':')
ax.set_yscale('log')
ax.set_ylim([2, 1e-8])
ax.set_xlabel('T [K]')
ax.set_ylabel('P [bar]')
ax.legend()
fig.savefig('/scratch/s2555875/plot/early_earth_TP_comp.png')
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
