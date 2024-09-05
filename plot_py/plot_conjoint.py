#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import pickle
import os
import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules
import parallel_functions as pf

nsim_dist = 15 # distances, for now
nsim_star = 13 # stars, for now

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output')
plot_folder = os.path.join(scratch, 'plot')
cfg_folder = os.path.join(output_folder, 'cfg')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))

min_flux_met = 1.058e-9
max_flux_met = 2.646e-8
#%%    
def get_surf_temp(file):
    surface_temperature = np.genfromtxt(file, dtype = None, skip_header=1, comments = '#', max_rows = 4, names = True)['Temp'][0]
    return surface_temperature

def check_convergence(dat):
    ''' It checks whether the convergence criteria has been met in the given simulation (dat is the 
        already read-in data). Template is taken from the VULCAN code (Tsai et al 2017, 2020).'''
    yconv_cri = 0.01
    yconv_min = 0.1
    slope_cri = 1.e-4
    longdy = dat['variable']['longdy']
    longdydt = dat['variable']['longdydt']
    slope_min = min( np.amin(dat['atm']['Kzz']/(0.1*dat['atm']['Hp'][:-1])**2) , 1.e-8)
    slope_min = max(slope_min, 1.e-10)
    if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min):
        return 'black' # True
    else:
        return 'none' # False
    
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
    ''' Calculates the rainout rate of the given species and returns it with units of kg/m2/yr.'''
    rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = rain_rate * 5.237e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

def plot_meshgrid(distance, teff, values, val_label, conv = None, figname = None, f_size = 13, l_size = 12,\
                  yscale = 'log', mol = 'HCN', norm = 'linear', met_flux = True):
    ''' Plots values, let it be rainout rate, end-of-simulation time, convergence, etc., in a pcolormesh plot.
        Diastance is X, teff is Y and values is C in documentation. Figure can be saved if needed.
        It is possible to highlight converged simulations with an edge using conv which shoukld be a 1D list of
        edgecolours.'''
    fig, ax = plt.subplots(tight_layout = True)
    a, T = np.meshgrid(distance, teff)
    cmap = plt.get_cmap('magma')
    cmap.set_under('none')
    cm = ax.pcolormesh(a, T, values, cmap = cmap, norm = norm, edgecolor = conv, linewidth = 0.05, snap = True)
    cbar = fig.colorbar(cm)
    cbar.set_label(val_label, fontsize = f_size)
    if met_flux:
        cbar.ax.axhline(min_flux_met, c = 'w')
        cbar.ax.axhline(max_flux_met, c = 'w')
    #for i in range(conv.shape[0]):
    #    for j in range(conv.shape[1]):
    #        if conv[i][j] == 'black':
    #            a_anchor = a[i][j]
    #            a_size = a[i, j + 1] - a[i, j]
    #            T_anchor = T[i][j]
    #            T_size = T[i + 1, j] - T[i, j]
    #            ax.add_patch(plt.Rectangle((a_anchor, T_anchor), a_size, T_size, fc='none', ec='black', lw=0.5, clip_on=False))
    ax.set_ylabel(r'T$_{eff}$ [K]', fontsize = f_size)
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]', fontsize = f_size)
    #ax.set_xscale('log')
    ax.invert_xaxis()
    ax.tick_params(which = 'both', direction = 'out', width = 1, length = 4)
    ax.tick_params(axis = 'both', labelsize = l_size)
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))


#%%
T_eff_list = list(star_df.T_eff)
a_list = []

for star in star_df.Name:
    #new_a_list = pf.semi_major_list_from_Seff(star_df, star, nsim_dist, factor = 1.1)
    new_a_list = pf.Seff_list(star_df, star, nsim_dist, factor = 1.1)
    i = 0
    for f in sorted(os.listdir(output_folder)):
        if 'star_' + star in f: # they are not in separate folder so dealing with only relevant files
            if new_a_list[i] not in a_list:
                a_list.append(new_a_list[i]) # save the new value
                i += 1
            
ncol = len(a_list)
#conv_matrix = [[0 for _ in range(ncol)] for _ in range(len(T_eff_list))] # matrix will be T_eff x a size (=ncol)
conv_matrix = [['none' for _ in range(ncol)] for _ in range(len(T_eff_list))] # matrix will be T_eff x a size (=ncol)
rain_matrix = [list(np.zeros(ncol)) for _ in range(len(T_eff_list))]
end_time_matrix = [list(np.zeros(ncol)) for _ in range(len(T_eff_list))]

for j,star in enumerate(star_df.Name):
    #new_a_list = pf.semi_major_list_from_Seff(star_df, star, nsim_dist, factor = 1.1)
    new_a_list = pf.Seff_list(star_df, star, nsim_dist, factor = 1.1)
    row_idx = j
    i = 0
    for f in sorted(os.listdir(output_folder)):
        if 'star_' + star in f: # they are not in separate folder so dealing with only relevant files
            column_idx = np.where(a_list == new_a_list[i])[0][0]
            with open(os.path.join(output_folder, f), 'rb') as handle:
                data = pickle.load(handle)
            
            conv_matrix[row_idx][column_idx] = check_convergence(data)
            rain_matrix[row_idx][column_idx] = rainout(data)
            end_time_matrix[row_idx][column_idx] = data['variable']['t']
            i += 1
#%%
#T_eff_list.append(T_eff_list[-1] + (T_eff_list[-1]-T_eff_list[-2]))
#a_max = max(a_list)
#a_max2 = sorted(a_list)[-2]
#a_list.append(a_max + (a_max-a_max2))
# %%
flat_conv = []
for row in conv_matrix:
    flat_conv.extend(row)
print('#converged: {}'.format(flat_conv.count('black')))
#%%
plot_meshgrid(a_list, T_eff_list, rain_matrix, r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', conv = flat_conv, norm = mc.LogNorm(vmin = 1e-20), figname = 'HCN_rainout_conjoint.pdf')
#%%
plot_meshgrid(a_list, T_eff_list, rain_matrix, r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', conv = flat_conv, norm = mc.LogNorm(vmin = 1e-20), figname = 'HCN_rainout_conjoint_S_eff.pdf')
# %%
plot_meshgrid(a_list, T_eff_list, end_time_matrix, 'End-of-simulation time [s]', conv = flat_conv, norm = mc.LogNorm(), figname = 'endtime_conjoint.pdf', met_flux = False)
# %%
plot_meshgrid(a_list, T_eff_list, end_time_matrix, 'End-of-simulation time [s]', conv = flat_conv, norm = mc.LogNorm(vmin = 1e-4), figname = 'endtime_conjoint_S_eff.pdf', met_flux = False)
# %%
