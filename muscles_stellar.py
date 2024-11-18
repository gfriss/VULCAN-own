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