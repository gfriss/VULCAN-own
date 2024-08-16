#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import parallel_functions as pf

nsim_dist = 15 # distances, for now
nsim_star = 13 # stars, for now

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output')
plot_folder = os.path.join(scratch, 'plot')
cfg_folder = os.path.join(output_folder, 'cfg')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))
#%%
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
    
def get_star_temp(star_name, df = star_df):
    return df.loc[df.Name == star_name].T_eff.iloc[0]

def get_star_name(file, df = star_df):
    star_name = ''
    if 'EARLY' in file:
        star_name = 'EARLY_SUN'
    else:
        star_name = file.split('_')[1]
    return star_name

def get_dist(file):
    cfg_file = os.path.join(cfg_folder, file[:-3] + 'txt') # from .vul to .txt
    orbit_r = 0.
    with open(cfg_file, 'r') as f:
        for line in f:
            if 'orbit_radius' in line:
                orbit_r = line.split()[2]
                break
    return orbit_r
    
def rainout(dat, rain_spec = 'HCN_rain', g_per_mol = 27):
    rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = rain_rate * 5.237e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

def plot_meshgrid(teff, distance, rain, figname = None, f_size = 13, l_size = 12, yscale = 'log', mol = 'HCN'):
    fig, ax = plt.subplots(tight_layout = True)
    # need a list of semi major axes, but it is different for each star so need to combine them ( and have corresponding T_eff list)
    # modify previous functions, make them into one and return the two lists...
    T, a= np.meshgrid(teff, distance)
    ax.pcolormesh(T, a, rain, cmap = 'magma')
    fig.colorbar()
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))


#%%
T_eff_list = list(star_df.T_eff)
a_list = []
non_zero_list = [] # list of lists

for i,f in enumerate(sorted(os.listdir(output_folder))):
    if 'star_' in f: # they are not in separate folder so dealing with only relevant files
        star = get_star_name(f)
        new_a = pf.semi_major_from_S_eff(star_df, star, factor = 1.1)
        
        if new_a not in a_list:
            a_list.append(new_a) # save the new value
            non_zero_list.append([i]) # make a new list for the current value
        else:
            a_idx = np.where(a_list == new_a)[0] # find the already existing value
            non_zero_list[a_idx].append(i) # extend the non-zero index for given value
        
#non_zero_list = [x for _,x in sorted(zip(a_list, non_zero_list))] # sorting the lists for tidiness
#a_list = sorted(a_list)
ncol = len(a_list)
conv_matrix = [list(np.zeros(ncol)) for _ in range(len(T_eff_list))] # matrix will be T_eff x a size
rain_matrix = [list(np.zeros(ncol)) for _ in range(len(T_eff_list))]

for i,f in enumerate(sorted(os.listdir(output_folder))):
    if 'star_' in f: # they are not in separate folder so dealing with only relevant files
        row_idx = i // nsim_dist # basically which star, should go until len(T_eff_list)
        column_idx = i % ncol # gives the current distance simulation number
        with open(os.path.join(output_folder, f), 'rb') as handle:
            data = pickle.load(handle)
        
        if i in non_zero_list[column_idx]:
            conv_matrix[] = check_convergence(data)
            rain_matrix[] = rainout(data)
        
# %%
