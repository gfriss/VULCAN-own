#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.lines import Line2D
import seaborn as sns
import pickle
import os
import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules
import parallel_functions as pf
from scipy.integrate import trapezoid
import pandas as pd
# setting up plot style
import plot_reset as pr
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6, grid = False)


nsim_dist = 15 # distances, for now
nsim_star = 13 # stars, for now
nsim_CtoO = 15 # C/O ratios, for now

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output')
plot_folder = os.path.join(scratch, 'plot')
cfg_folder = os.path.join(output_folder, 'cfg')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))
#network = ''
network = '_ncho'
nowash = '_nowash' # no washout case
archean_colour = 'k'
archean_file = os.path.join(scratch, 'output', 'archean'+network+nowash+'.vul')
with open(archean_file, 'rb') as handle:
    data_archean = pickle.load(handle)
    
min_mass_del = 4e20 * 1e-6 # g/yr, minimum mass delivery rate from kg/Gyr to g/yr
max_mass_del = 1e22 * 1e-6 # g/yr, maximum mass delivery rate from kg/Gyr to g/yr
met_hcn = 2472 # HCN meteoritic content in nmol/g by Smith et al. 2019
M_hcn = 27e-12 # kg/nmol, molar mass of HCN converted from 27 g/mol
A_earth = 5.1e14 # m2, surface area of the Earth
# meteoritic flux in kg/m2/yr
min_flux_met = met_hcn * min_mass_del * M_hcn / A_earth 
max_flux_met = met_hcn * max_mass_del * M_hcn / A_earth

#star_marker = {'M': 4, 'K': 5, 'G': 6, 'F': 7}
star_marker = {'M': 'o', 'K': '^', 'G': 'v', 'F': 's'}
star_markersize = [i*8+8 for i in range(6)] + [24, 32, 40] + [24, 32, 40] + [32] # separate for the M, K, G and F stars
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
    
def rainout(dat, rain_spec = 'HCN_rain', g_per_mol = 27):
    ''' Calculates the rainout rate of the given species and returns it with units of kg/m2/yr.'''
    #rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = trapezoid(dat['variable']['y_rain'][rain_spec], dat['atm']['zmco']) / dat['variable']['dt'] # 1/cm2s
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
        #cbar.ax.axhline(min_flux_met, c = 'w', lw = 2)
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
    cmap.set_over('none')
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

def get_seff(Temp, dist):
    ''' Calculates the effective stellar flux from the distance and temperature.'''
    llog = star_df.at[star_df.index[star_df['T_eff'] == Temp][0], 'L_log']
    seff = np.power(10, llog) / np.power(dist, 2)
    return seff

def get_stellar_type(Temp):
    ''' Gives back the stellar type.'''
    if Temp <= 4014: # values from table
        return 'M'
    elif Temp <= 5212:
        return 'K'
    elif Temp <= 6084:
        return 'G'
    else:
        return 'F'

def get_marker_size(Temp):
    ''' Calculates the marker size for the given effective stellar temperature.'''
    star_idx = np.where(star_df['T_eff'] == Temp)[0][0]
    return star_markersize[star_idx]

def get_marker_colour(idx, colours):
    ''' Calculates the marker colour for the given distances, but because the distances are different, it tracks the indexes.'''
    return colours[int((idx//nsim_CtoO)%nsim_dist)]
    
def read_pandas(file):
    dat = pd.read_csv(file, sep = '\t', skiprows=2, names = ['Teff', 'dist', 'CtoO', 'HCN_rain', 'H2O_rain'])
    # first sorting by Teff, dist and CtoO and dropping potential duplicates
    dat = dat.sort_values(by = ['Teff', 'dist', 'CtoO']).drop_duplicates(subset = ['Teff', 'dist', 'CtoO']).reset_index(drop = True)
    # then converting distance into effective stellar flux and replacing it
    dat['Seff'] = dat.apply(lambda row: get_seff(row['Teff'], row['dist']), axis = 1)
    # then calculating the marker tpye, size and colour
    dat['stellar_type'] = dat.apply(lambda row: get_stellar_type(row['Teff']), axis = 1)
    dat['marker_size'] = dat.apply(lambda row: get_marker_size(row['Teff']), axis = 1)
    colours = np.linspace(1, 0, nsim_dist, endpoint = True) # starts from the closest distance so brighetest colour
    dat.reset_index(inplace = True)
    dat['marker_colour'] = dat.apply(lambda row: get_marker_colour(row['index'], colours), axis = 1)
    return dat
    
def plot_3d(x, y, z, v, x_label, y_label, z_label, figsave = False):
    vmin = np.min(v[v>0])
    fig, ax = plt.subplots(tight_layout = True, subplot_kw={'projection': '3d'})
    cm = ax.scatter(x, y, z, c = np.array(v), norm = mc.LogNorm(vmin=vmin))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.invert_yaxis()
    cbar = fig.colorbar(cm)
    cbar.set_label(r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]')
    ax.view_init(elev=20, azim=-60, roll=0)
    
    if figsave:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/rain_3d'+network+'.pdf'), bbox_inches='tight')

def get_meshgrid_base(dat):
    base_CtoO = np.min(dat.CtoO) # the first C/O value in the dataframe
    x = np.array(dat.Seff[dat.CtoO == base_CtoO]).reshape(13,15)
    y = np.array(dat.Teff[dat.CtoO == base_CtoO]).reshape(13,15)
    z = np.array(dat.HCN_rain[dat.CtoO == base_CtoO]).reshape(13,15)
    e = np.full(shape=(13,15), fill_value='none')
    e[9,4] = archean_colour # the Archean Earth is always at the same position in the meshgrid
    e = e.flatten()
    return x, y, z, e

def get_meshgrid_many(dat, idx):
    ''' Takes the data from the pandas dataframe and returns the meshgrid for the base/minimum C/O value, as pandas dataframe loops through that first).
        It returns the Seff, Teff, HCN rainout rate and edge colours for the meshgrid plot.'''
    x = np.array(dat.Seff[dat.CtoO == dat.CtoO[idx]]).reshape(13,15)
    y = np.array(dat.Teff[dat.CtoO == dat.CtoO[idx]]).reshape(13,15)
    z = np.array(dat.HCN_rain[dat.CtoO == dat.CtoO[idx]]).reshape(13,15)
    e = np.full(shape=(13,15), fill_value='none')
    e[9,4] = archean_colour # the Archean Earth is always at the same position in the meshgrid
    e = e.flatten()
    return x, y, z, e
    
def plot_meshgrid_many(dat, val_label, figsave, rain_type = 'HCN_rain', met_flux = True):
    ''' Takes the data from the pandas dataframe, creates meshgrids for each C/O value and plots them in a pcolormesh plot.'''
    fig, axs = plt.subplots(ncols = 3, nrows = 5, tight_layout = True, layout = 'constrained', sharex = True, sharey = True)
    #fig.subplots_adjust(hspace = 0.3, wspace = 0.3)
    ax = axs.flatten()
    cmap = plt.get_cmap()
    cmap.set_under('none')
    vmin = np.nanmin(np.array(dat['HCN_rain'][dat['HCN_rain']>0])) # min value for the colour bar
    vmax = np.nanmax(np.array(dat['HCN_rain'])) # max value for the colour bar
    if met_flux:
        vmin = 0.9 * np.min([vmin, min_flux_met])
    ctoo = dat.CtoO.unique()
    for i,co in enumerate(ctoo):
        x, y, z, e = get_meshgrid_many(dat, i)
        cm = ax[i].pcolormesh(x, y, z, cmap = cmap, norm = mc.LogNorm(vmin=vmin, vmax=vmax))
        if i == 0: # the Archean simulation is for the first C/O value
            ax[i].pcolormesh(x, y, z, facecolors = 'none', edgecolors = e, lw = 2)
        if i%3 == 0: # y axis only for the first column
            ax[i].set_ylabel(r'T$_{eff}$ [K]')
        if i >= 12: # x axis only for the last row
            ax[i].set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
        ax[i].set_title('C/O = {}'.format(co))
        ax[i].invert_xaxis()
    cbar = fig.colorbar(cm, ax = axs, shrink = 0.8, location = 'right')
    cbar.set_label(val_label)
    if met_flux:
        cbar.ax.axhline(max_flux_met, c = 'grey', lw = 10)
    if figsave != None:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/'+rain_type+'out_many'+network+'.pdf'), bbox_inches='tight')
    
def plot_rainout_condensed(dat, figsave, rain_type = 'HCN_rain', met_flux = True):
    ''' Takes data from a pandas dataframe, then creates one plot that has the rainout rate information condensed.
        The y axis is the rainout rate (optionally with horizontal lines for the meteoritic delivery flux).
        The x axis is one of the three parameters (C/O), while marker size and colur are controlled by the other two parameters:
        The smaller the stellar effective temperature, the smaller the markersize + different marker for different stellar type; 
        the closer in the planet is, the brighter the colour.'''
    fig, ax = plt.subplots(tight_layout = True)
    seff_min, seff_max = np.min(dat.Seff), np.max(dat.Seff)
    plot_leg = [] 
    for st,marker in star_marker.items(): # loop through the stellar types and corresponding markers
        plot_leg.append(Line2D([0], [0], marker = marker, color = 'k', label = '{} star'.format(st), markersize = 6, linestyle = 'None'))
        dat_st = dat[dat['stellar_type'] == st]
        x = np.array(dat_st.CtoO)
        y = np.array(dat_st[rain_type])
        #ax.scatter(x, y, s = dat_st.marker_size, c = dat_st.marker_colour, marker = marker)
        cm = ax.scatter(x, y, s = dat_st.marker_size, c = dat_st.Seff, marker = marker, vmin = seff_min, vmax = seff_max)
    ax.set_ylabel(rain_type[:-5] + r' rainout [kg m$^{-2}$ yr$^{-1}$]')
    ax.set_xlabel('C/O')
    ax.set_yscale('log')
    ax.legend(handles = plot_leg, loc = 'upper right')
    cbar = fig.colorbar(cm)
    cbar.set_label(u'S$_{eff}$ [S$_\u2295$]')
    if met_flux:
        #ax.axhline(min_flux_met, c = 'k', lw = 2, ls = '--')
        ax.axhline(max_flux_met, c = 'gray', lw = 2, ls = '--')
    if figsave != None:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/'+rain_type+'_condensed'+network+'.pdf'), bbox_inches='tight')

def plot_hist(df, bins = 50, figsave = False, rain_type = 'HCN_rain', met_flux = True):
    ''' Takes the data from the pandas dataframe and plots a histogram of the rainout rates.'''
    #r = df[rain_type][df[rain_type] > 0] # only positive values
    #h, b = np.histogram(r, bins = bins)
    #log_b = np.logspace(np.log10(b[0]),np.log10(b[-1]),len(b))
    fig, ax = plt.subplots(tight_layout = True)
    #ax.hist(df[rain_type], bins = log_b)
    r = df[rain_type][df[rain_type] > 0] # only positive values
    r = r.dropna() # remove nan
    sns.histplot(r, bins = bins, stat = 'count', log_scale = True)
    #ax.set_xscale('log')
    ax.set_xlabel(rain_type[:-5] + r' rainout [kg m$^{-2}$ yr$^{-1}$]')
    ax.set_ylabel('Number of simulations')
    if met_flux:
        ax.axvline(max_flux_met, c = 'k', lw = 2, ls = '--', label = 'Maximum meteoritic flux')
    ax.axvline(1.835e-5, c = 'r', ls = '--', lw = 2, label = 'Archean Earth')
    ax.axvline(np.median(r), c = 'pink', ls = '--', lw = 2, label = 'Median')
    ax.legend()
    if figsave != None:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/'+rain_type+'_hist.pdf'), bbox_inches='tight')     

def plot_meshgrid_condensed(data, rain_type = 'HCN_rain', figsave = False):
    fig, ax = plt.subplots(tight_layout = True, ncols = 2)
    ax = ax.flatten()
    cmap = plt.get_cmap()
    cmap.set_under('none')
    # first doing the meshgrid for the base C/O value
    x, y, z, e = get_meshgrid_base(data)
    vmin = np.nanmin(z)
    vmin = 0.9 * np.min([vmin, min_flux_met])
    cm = ax[0].pcolormesh(x, y, z, cmap = cmap, norm = mc.LogNorm(vmin=vmin))
    ax[0].pcolormesh(x, y, z, edgecolors=e, facecolors='none', lw = 2)
    #ax[0].add_patch(plt.Rectangle(((x[9,4]+x[9,3])/2, y[9,4]), x[9,5]-x[9,4], y[10,4]-y[9,4], fc='none', edgecolor=archean_colour, lw=2, clip_on = False)) # Archean Earth
    cbar = plt.colorbar(cm, ax = ax[0], label = rain_type[:-5] + r' rainout [kg m$^{-2}$ yr$^{-1}$]')
    #cbar.set_label('{} rainout rate [kg m$^{-2}$ yr$^{-1}]'.format(rain_type[:-5]))
    ax[0].set_ylabel(r'T$_{eff}$ [K]')
    ax[0].set_xlabel(u'S$_{eff}$ [S$_\u2295$]')
    ax[0].invert_xaxis()
    # then the condensed plot for all C/O values
    seff_min, seff_max = np.min(data.Seff), np.max(data.Seff)
    plot_leg = [] 
    for st,marker in star_marker.items(): # loop through the stellar types and corresponding markers
        plot_leg.append(Line2D([0], [0], marker = marker, color = 'k', label = '{} star'.format(st), markersize = 6, linestyle = 'None'))
        dat_st = data[data['stellar_type'] == st]
        x = np.array(dat_st.CtoO)
        y = np.array(dat_st[rain_type])
        cm = ax[1].scatter(x, y, s = dat_st.marker_size, c = dat_st.Seff, marker = marker, vmin = seff_min, vmax = seff_max)
    ax[1].set_ylabel(rain_type[:-5] + r' rainout [kg m$^{-2}$ yr$^{-1}$]')
    ax[1].set_xlabel('C/O')
    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-12, None)
    ax[1].legend(handles = plot_leg, loc = 'upper right')
    plt.colorbar(cm, ax = ax[1], label = u'S$_{eff}$ [S$_\u2295$]')
    # add the meteoritic flux lines for both plots
    cbar.ax.axhline(max_flux_met, c = 'gray', lw = 2)
    ax[1].axhline(max_flux_met, c = 'gray', lw = 2, ls = '--', label = 'Maximum meteoritic flux')
    # add median and archean values for the second plot
    ax[1].axhline(np.nanmedian(z), c = 'pink', ls = '--', lw = 2, label = 'Median')
    ax[1].axhline(rainout(data_archean), c = archean_colour, lw = 2, ls = '--', label = 'Archean Earth')
    # and put letters on plots
    ax[0].text(0.03, 0.93, 'a)', transform=ax[0].transAxes)
    ax[1].text(0.03, 0.93, 'b)', transform=ax[1].transAxes)
    if figsave:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/'+rain_type+'meshgrid_condensed'+network+'.pdf'), bbox_inches='tight')

def plot_box_star(dat, rain_type = 'HCN_rain', figsave = False):
    type_stars = ['M', 'K', 'G', 'F']
    rains = []
    for t in type_stars:
        dt = dat[dat.stellar_type == t].dropna()
        rain_t = dt[rain_type][dt[rain_type] > 0]
        rains.append(rain_t)
        print('{} star planets:'.format(t))
        print('Average: {:.2e}\nMedian: {:.2e}\nStandard deviation: {:.2e}'.format(np.nanmean(rain_t), np.nanmedian(rain_t), np.nanstd(rain_t)))
    fig, ax = plt.subplots(tight_layout = True)
    bp = ax.boxplot(rains, labels = type_stars, patch_artist = True, notch = True)
    ax.set_yscale('log')
    ax.set_ylabel(rain_type[:-5] + r' rainout [kg m$^{-2}$ yr$^{-1}$]')
    ax.set_xlabel('Stellar type')
    for patch, c in zip(bp['boxes'], ['mistyrose', 'bisque', 'lightyellow', 'lightcyan']):
        patch.set_facecolor(c)
    if figsave:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/'+rain_type+'_stellar_types_boxplot'+network+'.pdf'), bbox_inches='tight')

def plot_box_CtoO(dat, rain_type = 'HCN_rain', figsave = False):
    type_CtoO = sorted(dat.CtoO.unique())
    rains = []
    for c in type_CtoO:
        dt = dat[dat.CtoO == c].dropna()
        rain_c = dt[rain_type][dt[rain_type] > 0]
        rains.append(rain_c)
        print('C/O = {}:'.format(c))
        print('Average: {:.2e}\nMedian: {:.2e}\nStandard deviation: {:.2e}'.format(np.nanmean(rain_c), np.nanmedian(rain_c), np.nanstd(rain_c)))
    fig, ax = plt.subplots(tight_layout = True)
    bp = ax.boxplot(rains, labels = type_CtoO, patch_artist = True, notch = True)
    ax.set_yscale('log')
    ax.set_ylabel(rain_type[:-5] + r' rainout [kg m$^{-2}$ yr$^{-1}$]')
    ax.set_xlabel('C/O')
    for patch, c in zip(bp['boxes'], sns.color_palette("magma_r", len(type_CtoO))):
        patch.set_facecolor(c)
    if figsave:
        fig.savefig(os.path.join(plot_folder,'rainout_rates/'+rain_type+'_CtoO_boxplot'+network+'.pdf'), bbox_inches='tight')
#%%
archean_rain = rainout(data_archean)
# 3D data
data_3d = read_pandas(os.path.join(scratch, 'star_dist_CtoO_rain.txt'))
#%% 3D plot
pr.reset_plt(ticksize = 11, fontsize = 13, fxsize = 9, fysize = 6, grid = False)
plot_3d(x = np.array(data_3d['Teff']), y = np.array(data_3d['Seff']), z = np.array(data_3d['CtoO']), v = np.array(data_3d['HCN_rain']), x_label = r'T$_{eff}$ [K]', y_label = u'S$_{eff}$ [S$_\u2295$]', z_label = 'C/O', figsave = True)
#%%
# meshdrid plots for base C/O value
Seff_mesh, Teff_mesh, HCN_rain_mesh, edge_mesh = get_meshgrid_base(data_3d)
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6, grid = False)
plot_meshgrid(Seff_mesh, Teff_mesh, HCN_rain_mesh, r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', edgec = sum(edge_mesh,[]), figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'.pdf')
plot_meshgrid_with_normed(Seff_mesh, Teff_mesh, HCN_rain_mesh, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed.pdf')
plot_meshgrid_with_normed(Seff_mesh, Teff_mesh, HCN_rain_mesh, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed_lognorm.pdf', norm = mc.LogNorm())
plot_meshgrid_with_normed(Seff_mesh, Teff_mesh, HCN_rain_mesh, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed_better.pdf', vmin = 1)
plot_meshgrid_with_normed(Seff_mesh, Teff_mesh, HCN_rain_mesh, archean_rain, 'Relative HCN rainout', figname = 'rainout_rates/HCN_rainout_conjoint_S_eff'+network+'_normed_worse.pdf', vmax = 1)
#%%
# meshgrid plot for all C/O values
pr.reset_plt(ticksize = 24, fontsize = 27, fxsize = 26, fysize = 30, grid = False)
plot_meshgrid_many(data_3d, val_label = r'HCN rainout [kg m$^{-2}$ yr$^{-1}$]', figsave = True, rain_type = 'HCN_rain')
#%%
# meshgrid and condensed plot side by side
pr.reset_plt(ticksize = 15, fontsize = 17, fxsize = 16, fysize = 6, grid = False)
plot_meshgrid_condensed(data_3d, rain_type = 'HCN_rain', figsave = True)
#%%
# condensed and histogram version
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6, grid = False)
plot_rainout_condensed(data_3d, figsave = True, rain_type = 'HCN_rain')
plot_hist(data_3d, figsave = True)
#%%
# probability plots (need normalisation!!! and flat prior...)
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6, grid = False)
p = gauss2d(Seff_mesh, Teff_mesh, 0.72, 5680, 0.072, 56.8) # figure out sigmas
plot_prior(Seff_mesh, Teff_mesh, p, 'Probability', edgec = sum(edge_mesh,[]), figname = 'probability/gaussian_prior.pdf')
plot_meshgrid_prob(Seff_mesh, Teff_mesh, HCN_rain_mesh, p, 'Probability', edgec = sum(edge_mesh,[]), figname = 'probability/gaussian_posterior.pdf')

p = cauchy2d(Seff_mesh, Teff_mesh, 0.72, 5680, 1) # figure out sigmas
plot_prior(Seff_mesh, Teff_mesh, p, 'Probability', edgec = sum(edge_mesh,[]), figname = 'probability/cauchy_prior.pdf')
plot_meshgrid_prob(Seff_mesh, Teff_mesh, HCN_rain_mesh, p, 'Probability', edgec = sum(edge_mesh,[]), figname = 'probability/cauchy_posterior.pdf')
#%%
# boxplots and statistics
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6, grid = False)
plot_box_star(data_3d, rain_type = 'HCN_rain', figsave = True)
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 12, fysize = 6, grid = False)
plot_box_CtoO(data_3d, rain_type = 'HCN_rain', figsave = True)
#%%