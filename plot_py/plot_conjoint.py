import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
#import parallel_functions as pf

nsim_dist = 15 # distances, for now
nsim_star = 13 # stars, for now

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output/star_dist/')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
#star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))

def check_hab_surf_temp(file):
    surface_temperature = np.genfromtxt(file, dtype = None, skip_header=1, comments = '#', max_rows = 4, names = True)['Temp'][0]
    if surface_temperature > 273. and surface_temperature < 373.:
        return True
    else:
        return False
    
def check_convergence(data):
    yconv_cri = 0.01
    yconv_min = 0.1
    slope_cri = 1.e-4
    longdy = data['variable']['longdy']
    longdydt = data['variable']['longdydt']
    slope_min = min( np.amin(data['atm']['Kzz']/(0.1*data['atm']['Hp'][:-1])**2) , 1.e-8)
    slope_min = max(slope_min, 1.e-10)
    if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min):
        return True
    else:
        return False
    
