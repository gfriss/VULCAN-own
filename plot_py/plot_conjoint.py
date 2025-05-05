#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import pickle
import os
import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules
import parallel_functions as pf
from scipy.integrate import trapezoid
# setting up plot style
import plot_reset as pr
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6, grid = False)


nsim_dist = 15 # distances, for now
nsim_star = 13 # stars, for now

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output')
plot_folder = os.path.join(scratch, 'plot')
cfg_folder = os.path.join(output_folder, 'cfg')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))
#network = ''
network = '_ncho'
archean_colour = 'k'
archean_file = os.path.join(scratch, 'output', 'archean'+network+'.vul')

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

def plot_meshgrid(x, y, values, val_label, edgec = 'none', figname = None, met_flux = True):
    ''' Plots values, let it be rainout rate, end-of-simulation time, convergence, etc., in a pcolormesh plot.
        Distance is X, Teff is Y. Figure can be saved if needed. Meteoritic flux can be added as well.'''
    vmin = np.min(values)
    if met_flux:
        vmin = 0.9 * np.min([vmin, min_flux_met])
    fig, ax = plt.subplots(tight_layout = True)
    cmap = plt.get_cmap()
    cmap.set_under('none')
    cm = ax.pcolormesh(x, y, values, cmap = cmap, norm = mc.LogNorm(vmin=vmin))
    ax.pcolormesh(x, y, values, facecolors='none', edgecolors=edgec, lw = 2)
    cbar = fig.colorbar(cm)
    cbar.set_label(val_label)
    if met_flux:
        cbar.ax.axhline(min_flux_met, c = 'w', lw = 2)
        cbar.ax.axhline(max_flux_met, c = 'w', lw = 2)
    ax.set_ylabel(r'T$_{eff}$ [K]')
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax.invert_xaxis()
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))
        
def plot_meshgrid_with_normed(x, y, values, norm_val, val_label, figname = None, vmin = None, vmax = None, norm = 'linear'):
    ''' Plots values that are normalised by the value of Archean Earth, let it be rainout rate, end-of-simulation time, convergence, etc., in a pcolormesh plot.
        Distance is X, Teff is Y. Figure can be saved if needed.'''
    fig, ax = plt.subplots(tight_layout = True)
    cmap = plt.get_cmap()
    cmap.set_under('none')
    cm = ax.pcolormesh(x, y, values/norm_val, cmap = cmap, vmin = vmin, vmax = vmax, norm = norm)
    cbar = fig.colorbar(cm)
    cbar.set_label(val_label)
    ax.set_ylabel(r'T$_{eff}$ [K]')
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax.invert_xaxis()
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))

def gauss2d(x, y, mux, muy, sigmax, sigmay):
    ''' 2D Gaussian function without correlation. x and y here are meshgrids like in the plot_meshgrid function.'''
    return (1/(2*np.pi*sigmax*sigmay)) * np.exp(-0.5 * ( ((x-mux)/sigmax)**2 + ((y-muy)/sigmay)**2 )) # this is the meshgrid version

def cauchy2d(x, y, x0, y0, gamma):
    ''' 2D Cauchy function. x and y here are meshgrids like in the plot_meshgrid function.'''
    return (gamma/(2*np.pi)) / ( (x-x0)**2 + (y-y0)**2 + gamma**2 )**1.5 # this is the meshgrid version

def normalise(x, y, values):
    ''' Normalises the values in the meshgrid using the trapz function from Scipy.
        First integration is in x, second is in y.'''
    Nx = [trapezoid(v[::-1], X[::-1]) for v,X in zip(values, x)] # Seff is flipped compared to orbital distance...
    N = trapezoid(Nx, y[:,0])
    return values/N
        
def plot_prior(x, y, values, val_label, edgec = 'none', figname = None):
    ''' Plots the prior used for the posterior.'''
    fig, ax = plt.subplots(tight_layout = True)
    cmap = plt.get_cmap()
    cmap.set_under('none')
    cm = ax.pcolormesh(x, y, normalise(x,y,values), cmap = cmap)
    ax.pcolormesh(x, y, values, facecolors='none', edgecolors=edgec, lw = 2)
    cbar = fig.colorbar(cm)
    cbar.set_label(val_label)
    ax.set_ylabel(r'T$_{eff}$ [K]')
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax.invert_xaxis()
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))
        
def plot_meshgrid_prob(x, y, values, prior, val_label, edgec = 'none', figname = None):
    ''' Plots values, let it be rainout rate, end-of-simulation time, convergence, etc., in a pcolormesh plot.
        Distance is X, Teff is Y. Figure can be saved if needed. Meteoritic flux can be added as well.'''
    fig, ax = plt.subplots(tight_layout = True)
    cmap = plt.get_cmap()
    cmap.set_under('none')
    cm = ax.pcolormesh(x, y, normalise(x,y,values*prior), cmap = cmap)
    ax.pcolormesh(x, y, values, facecolors='none', edgecolors=edgec, lw = 2)
    cbar = fig.colorbar(cm)
    cbar.set_label(val_label)
    ax.set_ylabel(r'T$_{eff}$ [K]')
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax.invert_xaxis()
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))

def plot_contour(x, y, values, val_label, figname = None, met_flux = True):
    ''' Plots values, let it be rainout rate, end-of-simulation time, convergence, etc., in a contour plot.
        Distance is X, Teff is Y. Figure can be saved if needed. Meteoritic flux can be added as well.'''
    vmin = np.min(values)
    if met_flux:
        vmin = 0.9 * np.min([vmin, min_flux_met])
    fig, ax = plt.subplots(tight_layout = True)
    cmap = plt.get_cmap()
    cmap.set_under('none')
    cm = ax.contourf(x, y, values, cmap = cmap, norm = mc.LogNorm(vmin=vmin))
    cbar = fig.colorbar(cm)
    cbar.set_label(val_label)
    if met_flux:
        cbar.ax.axhline(min_flux_met, c = 'w', lw = 2)
        cbar.ax.axhline(max_flux_met, c = 'w', lw = 2)
    ax.plot(0.72, 5680, color = archean_colour, marker = '*', markersize = 10)
    ax.set_ylabel(r'T$_{eff}$ [K]')
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax.invert_xaxis()
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))

def plot_tricontour(x, y, z, x_label, y_label, z_label, figname = None, met_flux = True):
    vmin = np.min(z)
    if met_flux:
        vmin = 0.9 * np.min([vmin, min_flux_met])
    fig, ax = plt.subplots(tight_layout = True)
    cmap = plt.get_cmap()
    cmap.set_under('none')
    cm = ax.tricontourf(x, y, z, cmap = cmap, norm = mc.LogNorm(vmin = vmin))
    cbar = fig.colorbar(cm)
    cbar.set_label(z_label)
    if met_flux:
        cbar.ax.axhline(min_flux_met, c = 'w', lw = 2)
        cbar.ax.axhline(max_flux_met, c = 'w', lw = 2)
    ax.plot(0.72, 5680, color = archean_colour, marker = '*', markersize = 10)
    ax.invert_xaxis()
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))
        
def plot_hist2d(x, y, weights, bins, weights_label, figname = None, met_flux = True):
    vmin = np.min(weights)
    if met_flux:
        vmin = 0.9 * np.min([vmin, min_flux_met])
    fig, ax = plt.subplots(tight_layout = True)
    cm = ax.hist2d(x=x, y=y, bins = bins, weights=weights, norm = mc.LogNorm(vmin = vmin))
    cbar = fig.colorbar(cm[3])
    cbar.set_label(weights_label)
    if met_flux:
        cbar.ax.axhline(min_flux_met, c = 'w', lw = 2)
        cbar.ax.axhline(max_flux_met, c = 'w', lw = 2)
    ax.plot(0.72, 5680, color = archean_colour, marker = '*') # try to highlight that square instead?
    ax.invert_xaxis()
    ax.set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax.set_ylabel(r'T$_{eff}$ [K]')
    if figname != None:
        fig.savefig(os.path.join(plot_folder, figname))
#%%
# meshgrid version
Teff_list, Seff_list, rain_matrix, end_time_matrix = [], [], [], []
edge_matrix = []
for star, a_min, a_max, Llog, T_eff in zip(star_df.Name, star_df.a_min, star_df.a_max, star_df.L_log, star_df.T_eff):
    distances = np.linspace(a_min, a_max, nsim_dist, endpoint=True)
    new_Seff_list = np.power(10, Llog) / np.power(distances, 2)
    Seff_list.append(new_Seff_list)
    Teff_list.append(np.ones(nsim_dist)*T_eff)
    rain_star, end_time_star = [], []
    edge_matrix.append(['none']*nsim_dist)
    if star == 'EARLY_SUN': # closest to the Archean Earth simulation
        edge_matrix[-1][4] = archean_colour
    # now find file and make up matrices
    for i in range(nsim_dist):
        if i < 10:
            sim = 'star_' + star + '_dist_0' + str(i) + network + '.vul'
        else:
            sim = 'star_' + star + '_dist_' + str(i) + network + '.vul'
        rain_rate = 0.
        end_time = 0.
        sim_file = os.path.join(output_folder, sim)
        with open(sim_file, 'rb') as handle:
            data = pickle.load(handle)
        rain_rate = rainout(data)
        end_time = data['variable']['t']
        if os.path.isfile(sim[:-4]+'_rerun.vul'):
            sim = sim[:-4]+'_rerun.vul'
            with open(os.path.join(output_folder, sim), 'rb') as handle:
                data = pickle.load(handle)
            rain_rate = rainout(data)
            end_time += data['variable']['t']
        rain_star.append(rain_rate)
        end_time_star.append(end_time)
    rain_matrix.append(np.array(rain_star))
    end_time_matrix.append(np.array(end_time_star))
Seff_list = np.array(Seff_list)
Teff_list = np.array(Teff_list)
rain_matrix = np.array(rain_matrix)
end_time_matrix = np.array(end_time_matrix)
#%%
plot_meshgrid(Seff_list, Teff_list, rain_matrix, r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', edgec = sum(edge_matrix,[]), figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'.pdf')
plot_meshgrid(Seff_list, Teff_list, end_time_matrix, 'End-of-simulation time [s]', edgec = sum(edge_matrix,[]), figname = 'end_time/endtime_conjoint_S_eff'+network+'.pdf', met_flux = False)
#%%
with open(archean_file, 'rb') as handle:
    data_archean = pickle.load(handle)
archean_rain = rainout(data_archean)
#%%
plot_meshgrid_with_normed(Seff_list, Teff_list, rain_matrix, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed.pdf')
plot_meshgrid_with_normed(Seff_list, Teff_list, rain_matrix, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed_lognorm.pdf', norm = mc.LogNorm())
plot_meshgrid_with_normed(Seff_list, Teff_list, rain_matrix, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed_better.pdf', vmin = 1)
plot_meshgrid_with_normed(Seff_list, Teff_list, rain_matrix, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed_worse.pdf', vmax = 1)
#%%
# probability plots (need normalisation!!! and flat prior...)
p = gauss2d(Seff_list, Teff_list, 0.72, 5680, 0.072, 56.8) # figure out sigmas
plot_prior(Seff_list, Teff_list, p, 'Probability', edgec = sum(edge_matrix,[]), figname = 'probability/gaussian_prior.pdf')
plot_meshgrid_prob(Seff_list, Teff_list, rain_matrix, p, 'Probability', edgec = sum(edge_matrix,[]), figname = 'probability/gaussian_posterior.pdf')

p = cauchy2d(Seff_list, Teff_list, 0.72, 5680, 1) # figure out sigmas
plot_prior(Seff_list, Teff_list, p, 'Probability', edgec = sum(edge_matrix,[]), figname = 'probability/cauchy_prior.pdf')
plot_meshgrid_prob(Seff_list, Teff_list, rain_matrix, p, 'Probability', edgec = sum(edge_matrix,[]), figname = 'probability/cauchy_posterior.pdf')

#%%
# contour version
# uses the same data as the meshgrid version
plot_contour(Seff_list, Teff_list, rain_matrix, r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_contour.pdf')
plot_contour(Seff_list, Teff_list, end_time_matrix, 'End-of-simulation time [s]', figname = 'end_time/endtime_conjoint_S_eff'+network+'_contour.pdf', met_flux = False)
#%%
# tricontour version
Seff_list = []
rain_list, end_time_list, Teff_list = [], [], []
for star,a_min,a_max,Llog,T_eff in zip(star_df.Name, star_df.a_min, star_df.a_max, star_df.L_log, star_df.T_eff):
    distances = np.linspace(a_min, a_max, nsim_dist, endpoint = True)
    new_Seff_list = np.power(10, Llog) / np.power(distances, 2)
    Seff_list += list(new_Seff_list)
    Teff_list += [T_eff]*nsim_dist
    for i in range(nsim_dist):
        if i < 10:
            sim = 'star_' + star + '_dist_0' + str(i) + network + '.vul'
        else:
            sim = 'star_' + star + '_dist_' + str(i) + network + '.vul'
        rain_rate = 0.
        end_time = 0.
        sim_file = os.path.join(output_folder, sim)
        with open(sim_file, 'rb') as handle:
            data = pickle.load(handle)
        rain_rate = rainout(data)
        end_time = data['variable']['t']
        if os.path.isfile(sim[:-4]+'_rerun.vul'):
            sim = sim[:-4]+'_rerun.vul'
            with open(os.path.join(output_folder, sim), 'rb') as handle:
                data = pickle.load(handle)
            rain_rate = rainout(data)
            end_time += data['variable']['t']
        rain_list.append(rain_rate)
        end_time_list.append(end_time)
#%%
plot_tricontour(Seff_list, Teff_list, rain_list, u'S$_{eff}$ [S$_\u2295$]', r'T$_{eff}$ [K]', r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_tricontour.pdf')        
plot_tricontour(Seff_list, Teff_list, end_time_list, u'S$_{eff}$ [S$_\u2295$]', r'T$_{eff}$ [K]', 'End-of-simulation time [s]', figname = 'end_time/endtime_conjoint_S_eff'+network+'_tricontour.pdf')        
# %%
# 2D histogram version
# uses the same data as the tricontour version
plot_hist2d(Seff_list, Teff_list, rain_list, 10, r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_hist2D.pdf')
plot_hist2d(Seff_list, Teff_list, end_time_list, 10, 'End-of-simulation time [s]', figname = 'end_time/endtime_conjoint_S_eff'+network+'_hist2D.pdf', met_flux = False)
# %%