'''
Script translating (already saved) MUSCLES spectra into stellar flux file acceptable for VULCAN.
Base code is taken from VULCAN (Tsai et al 2017.,2022.) directory.
The code scale the spectra back to surface of star.
'''

#%%
from astropy.io import fits
from astropy.constants import R_sun,pc
from astropy import units as u
import pandas as pd
import os
from collections import defaultdict
#%%
flux_dir = '/scratch/s2555875/stellar_flux'
fits_dir = os.path.join(flux_dir, 'fits_files')

r_sun = R_sun / u.m

types = defaultdict(lambda: 'float', Name = 'str', Type = 'str', T_eff_source = 'str', L_R_source = 'str', Dist_source = 'str')
stellar_table = pd.read_csv(os.path.join(flux_dir, 'stellar_params.csv'), dtype = types)
#%%
for filename in os.listdir(fits_dir):
    star = filename[:-5] # without .fits so only the name of the given star
    f = os.path.join(fits_dir, filename)
    spec = fits.getdata(f, memmap = True)
    new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'
    r_star = r_sun * stellar_table.loc[stellar_table.Name == star.upper()].R.iloc[0]
    d_star = stellar_table.loc[stellar_table.Name == star.upper()].Distance.iloc[0] * pc / u.m

    for wl,flux in zip(spec['WAVELENGTH'],spec['FLUX']):
        #if flux > 0: # there's negative values..., could do something smoother because it's still a big drop many a times
        new_str += '{:<12}'.format(wl*0.1) + "{:>12.2E}".format(float(flux*10. *(d_star/r_star)**2        )) + '\n'
        #else:
        #    new_str += '{:<12}'.format(wl*0.1) + "{:>12.2E}".format(0.5) + '\n' # trying whether zeros affect the convergence
    
    new_file = os.path.join(flux_dir,star + '.txt')
    with open(new_file, 'w') as g:
        g.write(new_str)
# %%
# testing spectra  to make sure which one VULCAN and HELIOS uses
# i.e. top of tmopsher or on surface of star...
import numpy as np
import matplotlib.pyplot as plt
# solar spectra from VULCAN:
gueymard_vulcan = np.genfromtxt('/home/s2555875/VULCAN-2/atm/stellar_flux/Gueymard_solar.txt', comments='#', names = ['lambda', 'flux'])
lam_vul = gueymard_vulcan['lambda']
F_vul = gueymard_vulcan['flux']
# solar spectra from Gueymard et al (2018) directly
# this is definitely at the top of the atmosphere
gueymard_og = np.genfromtxt('/scratch/s2555875/stellar_flux/gueymard_og.txt', comments='#', usecols=(0,1), names = ['lambda', 'flux'])
lam_og = gueymard_og['lambda']
F_og = gueymard_og['flux']*1e3 # unit conversion from Wm-2nm-1 to ergs-1cm-2nm-1
# plot for comparison with scaling original to solar surface
F_scaled = F_og*(1.496e8/6.955e5)**2 # F_TOA * (a/R_*)**2
plt.plot(lam_og, F_og, label = 'Gueymard et al. (2018)')
plt.plot(lam_vul, F_vul, label = 'from VULCAN')
plt.plot(lam_og, F_scaled, linestyle = '--', label = 'converted')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig('/scratch/s2555875/plot/solar_spectra_check.pdf')
#%%
# creating the early sun stellar spectrum from the file from Claire et al (2012) with a solar age of 0.6 Gyr

lam_f = np.genfromtxt('/scratch/s2555875/stellar_flux/early_sun_to_convert.txt', comments = '#', names = ['lambda', 'flux'])
# first convert to aangstrom because the conversiion factor uses aangstroms
lam = lam_f['lambda']*10
F = lam_f['flux']*0.1
# then do the conversion
F *= 1.98648e-8 / lam
# then convert back to nm
lam *= 0.1
F *= 10
# lastly, need to scale from top of atmosphere to solar surface
F *= (1.496e8/6.955e5)**2
# and write to file
new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'
for wl,flux in zip(lam,F):
    new_str += '{:<12}'.format(wl) + "{:>12.2E}".format(flux) + '\n'

new_file = os.path.join(flux_dir, 'early_sun.txt')
with open(new_file, 'w') as g:
    g.write(new_str)
#%%
import h5py
h5_dir = '/scratch/s2555875/HELIOS/star_tool/output'
flux_dir = '/scratch/s2555875/stellar_flux'
wl_conv = 1e7 # converting from cm to nm
flux_conv = 1e-7 # converting from erg/cm**3/s to erg/cm**2/s/nm
for filename in os.listdir(h5_dir):
    if 'sim' not in filename:
        continue
    sim = filename[:-3] # without .h5
    new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'
    f = h5py.File(os.path.join(h5_dir, filename), 'r')
    for wl, flux in zip(f['original']['phoenix']['lambda'][:], f['original']['phoenix'][sim][:]):
        new_str += '{:<12}'.format(wl*wl_conv) + "{:>12.2E}".format(flux*flux_conv) + '\n'
    f.close()
    with open(os.path.join(flux_dir, sim + '.txt'), 'w') as g:
        g.write(new_str)
# %%
