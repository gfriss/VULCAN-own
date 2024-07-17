#%%
import numpy as np
import sys

og_rad_file = 'atm/stellar_flux/Gueymard_solar.txt'
rad = np.genfromtxt(og_rad_file, comments = '#', skip_header = 1, dtype = None)
factor_A = [500, 60, 30, 20]
factor_B = [60, 10, 9, 7]
factor_def = 0.728 # fainter Sun by default, should affect the flux
if sys.argv[1] == 'A':
    factor = factor_A
elif sys.argv[1] == 'B':
    factor = factor_B
# %%
with open('atm/stellar_flux/Pearce_' + sys.argv[1] + '_solar_faint.txt', 'w') as f:
    f.write('# WL(nm)	 Flux(ergs/cm**2/s/nm)\n')
    for lam_flux in rad:
        if lam_flux[0] >= 0.1 and lam_flux[0] < 2.:
            f.write('{}\t{:.3e}\n'.format(lam_flux[0],lam_flux[1]*factor[0]))
        elif lam_flux[0] >= 2. and lam_flux[0] < 36.:
            f.write('{}\t{:.3e}\n'.format(lam_flux[0],lam_flux[1]*factor[1]))
        elif lam_flux[0] >= 36. and lam_flux[0] < 92.:
            f.write('{}\t{:.3e}\n'.format(lam_flux[0],lam_flux[1]*factor[2]))
        elif lam_flux[0] >= 92. and lam_flux[0] < 120.:
            f.write('{}\t{:.3e}\n'.format(lam_flux[0],lam_flux[1]*factor[3]))
        else:
            f.write('{}\t{:.3e}\n'.format(lam_flux[0],lam_flux[1]*factor_def))
# %%
