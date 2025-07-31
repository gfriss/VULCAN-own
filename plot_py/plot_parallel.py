#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import pickle
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
import os
import pandas as pd
wd = os.getcwd()
os.chdir('../')
from vulcan_cfg import yconv_cri, yconv_min, slope_cri, nl_ignore
os.chdir(wd)
# setting up plot style
import plot_reset as pr
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6)
# general parameters
nsim = 15 # hardcoded for now, change later...
out_folder = '/scratch/s2555875/output/'
plot_folder = '/scratch/s2555875/plot/'
sim_names = ['sim_0{}'.format(i) for i in range(nsim) if i < 10] + ['sim_{}'.format(i) for i in range(nsim) if i >= 10]
# end string of the VULCAN output files
end_str = {'BC': '_meteor', 'CtoO': '_CtoO', 'star': '_star', 'dist': '_dist'}
nowash = '_nowash' # no washout case
# setting chemical network for naming (empty if crahcno/no ignore as these are the defaults)
network = '_ncho'
# setting up the boundary condition case
bomb_rate = np.linspace(3e23, 1e25, nsim) # values from Pearce et al. (2022) Fig A4 min and max
bc_folder= '/scratch/s2555875/BC_files/'
bc_spec = 'H2'
bc_linestyle = '-'

# setting up the C/O case
co2_for_CtoO_range = np.linspace(0.1, 0.001, nsim, endpoint = True)
mixing_folder = '/scratch/s2555875/mixing_files/'

# setting up star case
star_df = pd.read_csv('/scratch/s2555875/stellar_flux/stellar_params.csv')
T_eff = star_df.T_eff
#T_eff = [2600+i*300 for i in range(nsim)]
#T_eff = [2800, 3100, 3400, 3700, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500]
# setting up the distance case
helios_output_folder = '/scratch/s2555875/HELIOS/output/'
stellar_spectra_folder = '/scratch/s2555875/stellar_flux/'

# setting up the distance case
a_list = np.linspace(0.839, 1.333, nsim, endpoint = True) #HZ limits from Kopprapau et al. (2013) are 0.99 and 1.7, let's explore a bit more, from Venus to 2 au

# base simulation of Archean
base_sim = out_folder+'archean'+network+nowash+'.vul'
with open(base_sim, 'rb') as handle:
    data_archean = pickle.load(handle)

# setting up plotting labels
archean_params = {'BC': 1.2e24 * 2.3/3.42, 'CtoO': 0.5143, 'star': 5680, 'dist': 1.}
xlab = {'BC': r'$\dot{M}_{del}$ [g/Gyr]', 'CtoO': 'C/O', 'star': r'T$_{eff}$ [K]', 'dist': 'Distance [AU]'}
xscale = {'BC': 'log', 'CtoO': 'linear', 'star': 'linear', 'dist': 'linear'}
legend_lab = {'BC': r'$\dot{{M}}_{{del}}$ = {:.2e} g/Gyr', 'CtoO': 'C/O = {:.3f}', 'star': r'T$_{{eff}}$ = {} K', 'dist': 'a = {:.3f} AU'}
archean_marker = 's' #'$\u2295$'
archean_colour = 'k'
# pressure levels for reaction rate plots
pressure_levels = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
# sim numbers to plot selected, detailed reaction rates, BC not needed
plot_idx = {'CtoO': ['archean', 4, 8, 12, 14], 'star': [0, 7, 'archean', 11, 12], 'dist': [0, 'archean', 8, 13, 14]}
plot_ls = ['-', '--', '-.', ':', (0, (1, 7))]
#%%
def lin(x, a, b):
    return a*x + b

def get_element_number(species, element):
    ''' Gets the number of an element in a given species. Used to calculate C/O ratio later.'''
    char_list = list(species)
    ele_number = 0
    for i,char in enumerate(char_list):
        if i < len(char_list)-1 and char == element:
            if char_list[i+1].isnumeric():
                ele_number += 1*int(char_list[i+1])
            else:
                ele_number += 1
        if i == len(char_list)-1 and char == element:
            ele_number += 1
    return ele_number

def calc_C_to_O(dat, mixing_file):
    ''' For simplicity and fewer calculations this funtion uses initial abundances, knowing the initial species that are none zero.
        It is generalised to get the species from the mixing ratio file used for the simulation.'''
    dat_species = dat['variable']['species']
    mix_data = np.genfromtxt(mixing_file, dtype = None, skip_header = 1, names = True, max_rows = 5) # max_rows so it wouldn't use too much unneccessary memory
    C_profile = np.zeros_like(dat['variable']['y_ini'][:,0])
    O_profile = np.zeros_like(dat['variable']['y_ini'][:,0])
    for name in mix_data.dtype.names:
        if name != 'Pressure' and 'C' in name:
            mul = get_element_number(name, 'C')
            C_profile += mul * dat['variable']['y_ini'][:, dat_species.index(name)]
        if name != 'Pressure' and 'O' in name:
            mul = get_element_number(name, 'O')
            O_profile += mul * dat['variable']['y_ini'][:, dat_species.index(name)]

    return np.sum(C_profile) / np.sum(O_profile)

def calc_mixing_h2(h2_bar):
    new_total = 1 + h2_bar
    return h2_bar / new_total

def read_in(sim_type, number_of_sim, start_str = ''):
    ''' Reads in all the data for a specific simulation type. The simulation type dictates how the files
        are found and what extra parameters this function returns. It builds on the idea that output
        files have a format of start_str+sim_+sim_number+sim_type+.vul.'''
    dat_list = [] # list to store the results (dicts)
    H2_flux = []
    T_surface = []
    ctoo = []
    for i in range(number_of_sim):
        sim = start_str + sim_names[i] # setting up the names of files
        sim += end_str[sim_type]
        sim_file = sim+'{}{}.vul'.format(network, nowash) # adding network and no washout if needed
        with open(out_folder+sim_file, 'rb') as handle:
            data_i = pickle.load(handle)
            dat_list.append(data_i)
        # if rerun was needed, combiune the results
        sim_rerun_file = os.path.join(out_folder,sim_file[:-4] + '_rerun.vul')
        if os.path.exists(sim_rerun_file):
            with open(sim_rerun_file, 'rb') as handle:
                data_rerun = pickle.load(handle)
            # change values that are carried over from first run, e.g. t, t_time, y_time
            data_rerun['variable']['t'] += data_i['variable']['t']
            data_rerun['variable']['t_time'] = np.concatenate((data_i['variable']['t_time'], data_rerun['variable']['t_time']+data_i['variable']['t']))
            data_rerun['variable']['y_time'] = np.concatenate((data_i['variable']['y_time'], data_rerun['variable']['y_time']), axis = 0)
            dat_list[i] = data_rerun
        # both networks use the same parameter file, so the name is the same
        # getting the parameters matching the run type
        if sim_type == 'BC':
            with open(bc_folder+'BC_bot_'+sim+'.txt') as f:
                for line in f:
                    lin = line.split()
                    if line[0] != '#' and lin[0] == bc_spec:
                        H2_flux.append(float(lin[1]))
                        break
        elif sim_type == 'CtoO':
            ctoo.append(calc_C_to_O(data_i, mixing_folder + sim + 'mixing.txt'))
        elif sim_type == 'dist':
            surface_temp = np.genfromtxt(helios_output_folder + sim + '/{}_tp.dat'.format(sim), dtype = None, skip_header = 2, usecols = (1))[0]
            T_surface.append(surface_temp)

    if sim_type == 'BC':
        return dat_list, H2_flux
    elif sim_type == 'CtoO':
        return dat_list, ctoo
    elif sim_type == 'dist':
        return dat_list, T_surface
    else:
        return dat_list

def rainout(dat, rain_spec = 'HCN_rain', g_per_mol = 27):
    ''' Calculates the rainout rate of the given species and returns it with units of kg/m2/yr.'''
    #rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = trapezoid(y=dat['variable']['y_rain'][rain_spec], x=dat['atm']['zmco']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = rain_rate * 5.237e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)

def plot_rain(hcn_rain_list, param_list, sim_type, figsave, extra_list = [], yscale = 'log', rain_spec = 'HCN_rain'):
    ''' Function to plot rainout rates of different simulatio types along with the Archean results.'''
    fig, ax = plt.subplots(tight_layout = True)
    ax.plot(param_list, hcn_rain_list, color = 'navy', linestyle = '', marker = '.', markersize = 10)
    ax.plot(archean_params[sim_type], rainout(data_archean, rain_spec = rain_spec), color = archean_colour, marker = archean_marker, markersize = 10)
    ax.set_xlabel(xlab[sim_type])
    ax.set_xscale(xscale[sim_type])
    ax.set_ylabel(rain_spec[:-5] + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]')
    ax.set_yscale(yscale)
    if sim_type == 'BC': # plotting comparison and H2 flux
        ax1 = ax.twiny()
        ax1.plot(extra_list, hcn_rain_list, linestyle = '')
        ax1.set_xlabel(r'H$_2$ flux [cm$^{-2}$ s$^{-1}$]')
        ax1.set_xscale('log')
    
    elif sim_type == 'dist':
        ax1 = ax.twiny()
        ax1.plot(extra_list, hcn_rain_list, linestyle = '')
        ax1.set_xlabel(r'T$_{surf}$ [K]')
        ax1.invert_xaxis()

    if figsave:
        fig.savefig(plot_folder + 'rainout_rates/'+rain_spec+'out'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def calc_rate(dat, spec_list, re_id, n):
    ''' Calculates the reaction rate by multiplying the reaction coefficient with the number densty of
        the reagents.'''
    rate = dat['variable']['k'][re_id][n]
    for sp in spec_list:
        if sp != 'M':
            rate *= dat['variable']['y'][n, dat['variable']['species'].index(sp)]
        else:
            rate *= dat['atm']['n_0'][n]
    return rate

def get_species(eq_side):
    ''' Returns the species in a given reaction in an array.'''
    side_split = eq_side.split('+')
    if len(side_split) == 1: # stripping them from white spaces
        side_split = np.array([side_split[0].strip()]) # with array length 1 the other method fails so doing it separately
    else:
        side_split = np.array([r.strip() for r in side_split])
    return side_split

def print_max_re(dat, mol, n = 0, prod = False):
    ''' Finds and prints the reaction with the highest reaction rate for a given molecule.
        This reaction can be production or destruction reaction, controlled by the prod parameter
        (True and False, respectively).'''
    reaction = ''
    k_rea = 0
    for k,v in dat['variable']['Rf'].items():
        reagents_spec, products_spec = [], []
        re_rate = 0
        reagents_products = v.split('->') # separating production and destruction
        reagents_spec = get_species(reagents_products[0])
        products_spec = get_species(reagents_products[1])
        if prod == False and mol in reagents_spec: # destruction so it is on the reagents side
            re_rate = calc_rate(dat, reagents_spec, k, n)
        elif prod == True and mol in products_spec: # production side
            re_rate = calc_rate(dat, products_spec, k, n)
        if re_rate > k_rea:
            k_rea = re_rate
            reaction = v
    print('The highest reaction rate in layer {} is\nk = {:.3e} cm-3s-1\nfor the reaction {}.'.format(n, k_rea, reaction))

def print_max_n_re(dat, mol, first_n, n = 0, prod = False):
    ''' Similar to previous function but it finds the first n reactions.'''
    reaction = np.empty(first_n, dtype = np.dtype('U42')) # extended the allowed character number to 42...
    reaction = np.array(reaction)
    k_rea = np.zeros(first_n)
    for k,v in dat['variable']['Rf'].items():
        re_rate = 0.
        reagents_spec, products_spec = [], []
        reagents_products = v.split('->') # separating production and destruction
        reagents_spec = get_species(reagents_products[0]) 
        products_spec = get_species(reagents_products[1])
        if (prod == False) and np.any(reagents_spec == mol): # destruction so it is on the reagents side
            re_rate = calc_rate(dat, reagents_spec, k, n)
        elif (prod == True) and np.any(products_spec == mol): # production side
            re_rate = calc_rate(dat, products_spec, k, n)
        if np.any(re_rate > k_rea):
            position = np.where(re_rate > k_rea)[0][0] # k_re should be decreasing, so need first intance when new reaction rate is greater
            k_roll = np.roll(k_rea[position:], 1) # shifting value to the right starting from position of new value
            k_roll[0] = re_rate # replace the new first one (previous lowest) with the new rate
            k_rea = np.concatenate((k_rea[:position], k_roll))
            reaction_roll = np.roll(reaction[position:], 1)
            reaction_roll[0] = v
            reaction = np.concatenate((reaction[:position], reaction_roll))
    print('The highest {} reaction rates in layer {} are k = {} cm-3s-1 for the reaction {}.'.format(first_n, n, k_rea, reaction))
    
def plot_vertical_n(dat_list, spec, param_list, sim_type, figsave):
    ''' Plots the vertical profile of the given species (spec) at the end of simulation for a list
        of simulations.'''
    fig, ax = plt.subplots(tight_layout = True)
    min_mixing = 1
    for d,p in zip(dat_list, param_list):
        ax.plot(d['variable']['ymix'][:, d['variable']['species'].index(spec)], d['atm']['pco']/1e6, label = legend_lab[sim_type].format(p))
        min_mixing = np.min([min_mixing, np.min(d['variable']['ymix'][:, d['variable']['species'].index(spec)])])  
    min_mixing = np.max([0.8*min_mixing, 1e-15])
    ax.set_xscale('log')
    ax.set_xlabel('Mixing ratio of {}'.format(spec))
    ax.set_xlim(min_mixing, None)
    ax.set_ylabel('Pressure [bar]')
    ax.set_yscale('log')
    ax.invert_yaxis()
    fig.legend(bbox_to_anchor = (1.35, 0.97))
    if figsave:
        fig.savefig(plot_folder + 'vertical_profiles/'+spec+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def plot_vertical_many(dat_list, param_list, sim_type, figsave, species_list = ['CH4', 'O', 'OH', 'CN', 'HNCO', 'H2CN', 'C2H3', 'C2H3CN', 'C2H6']):
    fig, ax = plt.subplots(tight_layout = True)
    ax.set_prop_cycle(color = plt.get_cmap('tab10').colors)
    plot_leg = [] # linestyle legends
    for i,idx in enumerate(plot_idx[sim_type]):
        if idx == 'archean': # if it is the archean simulation, use the archean data
            plot_leg.append(Line2D([], [], linestyle = plot_ls[i], color = 'black', label = legend_lab[sim_type].format(archean_params[sim_type]) + ' (Archean)')) # Archean legend
        else:
            plot_leg.append(Line2D([], [], linestyle = plot_ls[i], color = 'black', label = legend_lab[sim_type].format(param_list[idx]))) # other legends
    for i,idx in enumerate(plot_idx[sim_type]):
        if idx == 'archean': # if it is the archean simulation, use the archean data
            d = data_archean
        else:
            d = dat_list[idx]
        for sp in species_list:
            p = ax.plot(d['variable']['ymix'][:, d['variable']['species'].index(sp)], d['atm']['pco']/1e6, linestyle = plot_ls[i])
            if i == 0: # only first loop needs to be added to the legend
                plot_leg.append(Line2D([], [], linestyle = '-', color = p[0].get_color(), label = sp)) # species legends
        ax.set_prop_cycle(color = plt.get_cmap('tab10').colors) # resetting the colour cycle
    ax.set_xscale('log')
    ax.set_xlabel('Mixing ratio of species')
    ax.set_xlim(1e-20, 1e-1)
    ax.set_ylabel('Pressure [bar]')
    ax.set_yscale('log')
    ax.invert_yaxis()
    fig.legend(handles = plot_leg, bbox_to_anchor = (1.34, 0.87))
    if figsave:
        fig.savefig(plot_folder + 'vertical_profiles/many'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')
        
def plot_end_time(dat_list, figsave, sim_type = None):
    ''' Plots the end-of-simulation times for a list of simulations. It is used as a way of 
        seeing what difference there are between both converged and not converged simulations.'''
    fig, ax = plt.subplots(tight_layout = True)
    for i in range(len(dat_list)):
        ax.plot(i, dat_list[i]['variable']['t'], linestyle = '', marker = 'o', color = 'red')
    ax.set_yscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Simulation number')
    ax.set_ylabel('End of simulation time [s]')
    if figsave:
        fig.savefig(plot_folder + 'end_time/'+end_str[sim_type][1:]+network+'.pdf', bbox_inches = 'tight')

def plot_evo_layer(dat_list, param_list, spec, layer, figsave, sim_type = None):
    ''' Plots the evolution of a given species in a given layer for a list of simulations.'''
    fig, ax = plt.subplots(tight_layout = True)
    for d,p in zip(dat_list, param_list):
        ax.plot(d['variable']['t_time'], d['variable']['y_time'][:, layer, d['variable']['species'].index(spec)], label = legend_lab[sim_type].format(p))
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('n [cm-3]')
    ax.set_yscale('log')
    ax.set_ylim(1e-2,None)
    fig.legend(bbox_to_anchor = (1.35, 0.97))
    if figsave:
        fig.savefig(plot_folder + 'evolution/'+spec+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def check_convergence(dat):
    ''' It checks whether the convergence criteria has been met in the given simulation (dat is the 
        already read-in data). Template is taken from the VULCAN code (Tsai et al 2017, 2020).'''
    longdy = dat['variable']['longdy']
    longdydt = dat['variable']['longdydt']
    slope_min = min( np.amin(dat['atm']['Kzz']/(0.1*dat['atm']['Hp'][:-1])**2) , 1.e-8)
    slope_min = max(slope_min, 1.e-10)
    if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min):
        return True
    else:
        return False

def plot_convergence(dat_list, figsave, sim_type = None):
    ''' It checks and plots whether the convergence criteria has been met in all simulations in the given list. 
        Template for calculation is taken from the VULCAN code (Tsai et al 2017, 2020).'''
    fig, ax = plt.subplots(tight_layout = True)
    for i,d in enumerate(dat_list):
        if check_convergence(d):
            ax.plot(i, 'Yes', 'ro')
        else:
            ax.plot(i, 'No', 'ro')
    ax.set_xlabel('Simulation number')
    ax.set_ylabel('Converged')
    
    if figsave:
        fig.savefig(plot_folder + 'convergence/'+end_str[sim_type][1:]+network+'.pdf', bbox_inches = 'tight')

def plot_rain_converged(dat_list, rain_list, param_list, sim_type, figsave, rain_spec = 'HCN_rain', extra_list = [], plot_non_conv = False):
    ''' Plots the rainout rates for a list of simulations for a given simulation types. It 
        distinguishes between converged and non-converged simulations (full and empty circles, respectively).
        Plotting and convergence caalculations are taken from previous functions.'''
    file_start = 'conv_' # to distinguish between converged and non-converged files
    conv_rain_list, conv_param_list, conv_extra_list = [], [], []
    non_conv_rain_list, non_conv_param_list, non_conv_extra_list = [], [], []
    for i,d in enumerate(dat_list):
        if check_convergence(d):
            conv_rain_list.append(rain_list[i])
            conv_param_list.append(param_list[i])
            if extra_list:
                conv_extra_list.append(extra_list[i])
        else:
            non_conv_rain_list.append(rain_list[i])
            non_conv_param_list.append(param_list[i])
            if extra_list:
                non_conv_extra_list.append(extra_list[i])
    gpm = {'HCN_rain': 27, 'H2O_rain': 18} # g per mol for the species
    rain_archean = rainout(data_archean, rain_spec = rain_spec, g_per_mol = gpm[rain_spec])
    
    popt, pcov = curve_fit(lin, conv_param_list, conv_rain_list)
    param_x = np.linspace(conv_param_list[0], conv_param_list[-1], 42, endpoint = True)
    fig, ax = plt.subplots(tight_layout = True)
    ax.plot(conv_param_list, conv_rain_list, color = 'navy', linestyle = '', marker = '.', markersize = 10)
    if plot_non_conv and non_conv_param_list: # told to plot non converged and there are non converged sims
        ax.plot(non_conv_param_list, non_conv_rain_list, color = 'navy', linestyle = '', marker = '.', fillstyle = 'none', markersize = 10)
        file_start += 'non_conv_'
    ax.plot(archean_params[sim_type], rain_archean, color = archean_colour, marker = archean_marker, markersize = 10)
    ax.plot(param_x, lin(param_x, popt[0], popt[1]), linestyle = '--', alpha = 0.4, color = 'r')
    ax.set_yscale('log')
    ax.set_xscale(xscale[sim_type])
    ax.set_xlabel(xlab[sim_type])
    if sim_type == 'BC':
        ax1 = ax.twiny()
        ax1.plot(conv_extra_list, conv_rain_list, linestyle = '')
        if plot_non_conv and non_conv_param_list: # told to plot non converged and there are non converged sims
            ax1.plot(non_conv_extra_list, non_conv_rain_list, color = 'navy', linestyle = '', marker = '.', markerfacecolor = 'none')
        ax1.set_xlabel(r'H$_2$ flux [cm$^{-2}$ s$^{-1}$]')
        ax1.set_xscale('log')
    elif sim_type == 'dist':
        ax1 = ax.twiny()
        ax1.plot(conv_extra_list, conv_rain_list, linestyle = '')
        if plot_non_conv and non_conv_param_list: # told to plot non converged and there are non converged sims
            ax1.plot(non_conv_extra_list, non_conv_rain_list, color = 'navy', linestyle = '', marker = '.', markerfacecolor = 'none')
        ax1.set_xlabel(r'T$_{surf}$ [K]')
        ax1.invert_xaxis()

    ax.set_ylabel(rain_spec[:-5] + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]')

    if figsave:
        fig.savefig(plot_folder + 'rainout_rates/'+file_start+rain_spec+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')
    
def get_species(eq_side):
    ''' Returns the species in a given reaction side in a list.'''
    side_split = eq_side.split('+')
    side_split = [r.strip() for r in side_split] # stripping them from white spaces
    return side_split
    
def get_total_reaction_rate(dat, diag_sp = 'HCN'):
    species = dat['variable']['species']
    total_rate = np.zeros_like(dat['atm']['pco'])
    prod_rate = np.zeros_like(dat['atm']['pco'])
    dest_rate = np.zeros_like(dat['atm']['pco'])
    for re_id,rea in dat['variable']['Rf'].items():
        reagents_products = rea.split('->')
        reagents = get_species(reagents_products[0])
        products = get_species(reagents_products[1])
        if diag_sp not in reagents and diag_sp not in products: # if the species is not in the reaction, skip it
            continue
        forward_rate = dat['variable']['k'][re_id].astype(float)
        reverse_rate = dat['variable']['k'][re_id+1].astype(float)
        for sp in reagents: # calculating forward rate
            if sp == 'M': forward_rate *= dat['atm']['n_0']
            else: forward_rate *= dat['variable']['y'][:,species.index(sp)]
        for sp in products: # calculating reverse rate
            if sp == 'M': reverse_rate *= dat['atm']['n_0']
            else: reverse_rate *= dat['variable']['y'][:,species.index(sp)]
        
        if diag_sp in reagents: 
            dest_rate += np.array(forward_rate) # it gets destroyed when reactiong with something
            prod_rate += np.array(reverse_rate) # it is produced in the reverse reaction
        elif diag_sp in products: 
            prod_rate += np.array(forward_rate) # it is produced in the forward reaction
            dest_rate += np.array(reverse_rate) # it is destroyed when reacting with something
    total_rate = prod_rate - dest_rate
    return prod_rate, dest_rate, total_rate    
    
def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)    
    
def plot_tot_rate(dat_list, param_list, sim_type, diag_sp, figsave):
    fig, ax = plt.subplots(tight_layout = True)
    for d,p in zip(dat_list, param_list):
        pressure = d['atm']['pco']/1e6
        _, _, tot_rate = get_total_reaction_rate(d, diag_sp)
        ax.plot(tot_rate, pressure, label = legend_lab[sim_type].format(p))
        
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_xscale('symlog')
    ax.set_ylabel('Pressure [bar]')
    ax.set_xlabel(r'k$_{tot}$ [cm$^3$s$^{-1}$]')  
    fig.legend(bbox_to_anchor = (1.35,0.97))
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/total'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')
        
def plot_tot_rate_selected(dat_list, param_list, sim_type, diag_sp, figsave):
    fig, ax = plt.subplots(tight_layout = True)
    ax.set_prop_cycle(color = plt.get_cmap('tab10').colors, linestyle = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']) 
    for i,idx in enumerate(plot_idx[sim_type]):
        if idx == 'archean': # if it is the archean simulation, use the archean data
            d = data_archean
            leg = legend_lab[sim_type].format(archean_params[sim_type]) + ' (Archean)' # adding Archean to the legend
        else:
            d = dat_list[idx]
            leg = legend_lab[sim_type].format(param_list[idx])
        pressure = d['atm']['pco']/1e6
        _, _, tot_rate = get_total_reaction_rate(d, diag_sp)
        ax.plot(tot_rate, pressure, linestyle = plot_ls[i], label = leg)
        
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_xscale('symlog')
    ax.set_ylabel('Pressure [bar]')
    ax.set_xlabel(r'k$_{tot}$ [cm$^3$s$^{-1}$]')  
    ax.legend(loc = 'lower right')
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/selected_total'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')
        
def plot_prod_dest(dat_list, param_list, sim_type, diag_sp, figsave):
    fig, ax = plt.subplots(tight_layout = True, nrows = 2, ncols = 1, sharex = True, figsize = (6,8))
    ax = ax.flatten()
    for d,p in zip(dat_list, param_list):
        pressure = d['atm']['pco']/1e6
        p_rate, d_rate, _ = get_total_reaction_rate(d, diag_sp)
        ax[0].plot(p_rate, pressure, label = legend_lab[sim_type].format(p))
        ax[1].plot(d_rate, pressure)
    for i in range(len(ax)):
        ax[i].invert_yaxis()
        ax[i].set_yscale('log')
        ax[i].set_xscale('log')
        ax[i].set_ylabel('Pressure [bar]')
    ax[0].set_xlabel(r'k$_{prod}$ [cm$^3$s$^{-1}$]')
    ax[1].set_xlabel(r'k$_{dest}$ [cm$^3$s$^{-1}$]')  
    fig.legend(bbox_to_anchor = (1.45,0.85))
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/net'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def plot_prod_dest_selected(dat_list, param_list, sim_type, diag_sp, figsave):
    fig, ax = plt.subplots(tight_layout = True)
    plot_leg = []
    for i,idx in enumerate(plot_idx[sim_type]):
        if idx == 'archean':
            plot_leg.append(Line2D([], [], linestyle = plot_ls[i], color = 'black', label = legend_lab[sim_type].format(archean_params[sim_type]) + ' (Archean)')) # Archean legend
        else:
            plot_leg.append(Line2D([], [], linestyle = plot_ls[i], color = 'black', label = legend_lab[sim_type].format(param_list[idx])))
    plot_leg += [Line2D([], [], linestyle = '-', color = 'red', label = 'Production'), Line2D([], [], linestyle = '-', color = 'black', label = 'Destruction')] # species legends
    for i,idx in enumerate(plot_idx[sim_type]):
        if idx == 'archean': # if it is the archean simulation, use the archean data
            d = data_archean
        else:
            d = dat_list[idx]
        pressure = d['atm']['pco']/1e6
        p_rate, d_rate, _ = get_total_reaction_rate(d, diag_sp)
        ax.plot(p_rate, pressure, color = 'r', linestyle = plot_ls[i])
        ax.plot(d_rate, pressure, color = 'k', linestyle = plot_ls[i])    
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Pressure [bar]')
    ax.set_xlabel(r'k [cm$^3$s$^{-1}$]')
    fig.legend(handles=plot_leg, bbox_to_anchor = (0.47,0.97))
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/selected_net'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')
        
def plot_prod_dest_layer(dat_list, param_list, sim_type, layer, diag_sp, figsave):
    fig, ax = plt.subplots(tight_layout = True)
    prod, dest = [], []
    for d in dat_list:
        p_rate, d_rate, _ = get_total_reaction_rate(d, diag_sp)
        prod.append(p_rate[layer])
        dest.append(d_rate[layer])
    ax.plot(param_list, prod, label = 'Production', c = 'r')
    ax.plot(param_list, dest, label = 'Destruction', c = 'k')
    ax.set_xlabel(xlab[sim_type])
    ax.set_xscale(xscale[sim_type])
    ax.set_ylabel(r'k [cm$^{-3}$s$^{-1}$]')    
    ax.set_yscale('log')
    ax.legend()
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/layer_'+str(layer)+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def get_layer(dat, pressure):
    ''' Returns the layer that is the closest to the given pressure (in bar) in the given VULCAN data.'''
    return np.argmin(abs(dat['atm']['pco']/1e6 - pressure))

def plot_prod_dest_many_layer(dat_list, param_list, sim_type, pressures, ncols, diag_sp, figsave):
    n = len(pressures)
    nrows = n//ncols # setting up number of rows given columns
    if n%ncols != 0: # if not divisible by ncols, add one row
        nrows += 1
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize = (5*ncols,4*nrows), sharex = True, sharey = True, tight_layout = True)
    ax = ax.flatten()
    prod, dest = [], []
    # first get the rates for all simulations (list, by pressures, of lists, by simulations)
    for p in pressures:
        prod_p, dest_p = [], []
        for d in dat_list:
            p_rate, d_rate, _ = get_total_reaction_rate(d, diag_sp)
            prod_p.append(p_rate[get_layer(d, p)])
            dest_p.append(d_rate[get_layer(d, p)])
        prod.append(prod_p)
        dest.append(dest_p)
    
    for i,axes in enumerate(ax):
        axes.plot(param_list, prod[i], c = 'r')
        axes.plot(param_list, dest[i], c = 'k')
        if i%ncols == 0: # only first column gets y label
            axes.set_ylabel(r'k [cm$^{-3}$s$^{-1}$]')
        if i >= (nrows-1)*ncols: # only last row gets x label
            axes.set_xlabel(xlab[sim_type])
        axes.set_xscale(xscale[sim_type])
        axes.set_yscale('log')
        axes.legend(title = 'P = {:.1e} bar'.format(pressures[i]), loc = 'upper left')
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/many_layers'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def get_prod_dest_rates(dat, diag_sp):
    ''' Function to get the production and destruction rates for a given species (diag_sp) in a given, already read-in VULCAN simulation (dat).
        It returns the rates and the important reactions that are >10% of the total rate in two dictionaries that both have sub-dictionaries
        for production and destruction. The rates subdictionaries also have the total rate and the min and max xlim values for plotting.'''
    species = dat['variable']['species']
    rates = {'Production': {'xlim_min': 1e-1, 'xlim_max': 1e3, 'total': 0.}, 'Destruction': {'xlim_min': 1e-1, 'xlim_max': 1e3, 'total': 0.}}
    important_rates = {'Production': [], 'Destruction': []}
    tot_prod_rate = np.zeros_like(dat['atm']['pco'])
    tot_dest_rate = np.zeros_like(dat['atm']['pco'])
    for re_id,rea in dat['variable']['Rf'].items():
        reagents_products = rea.split('->')
        reagents = get_species(reagents_products[0])
        products = get_species(reagents_products[1])
        if diag_sp not in reagents and diag_sp not in products: # if the species is not in the reaction, skip it
            continue
        forward_rate = dat['variable']['k'][re_id].astype(float)
        reverse_rate = dat['variable']['k'][re_id+1].astype(float)
        rev_rea = ' -> '.join((' + '.join(products), ' + '.join(reagents))) # reverse reaction to be used as key in rates
        for sp in reagents: # calculating forward rate
            if sp == 'M': forward_rate *= dat['atm']['n_0']
            else: forward_rate *= dat['variable']['y'][:,species.index(sp)]
        for sp in products: # calculating reverse rate
            if sp == 'M': reverse_rate *= dat['atm']['n_0']
            else: reverse_rate *= dat['variable']['y'][:,species.index(sp)]
        
        if diag_sp in reagents:
            rates['Destruction'][rea] = forward_rate # it gets destroyed when reactiong with something
            tot_dest_rate += np.array(forward_rate)
            rates['Production'][rev_rea] = reverse_rate # it is produced in the reverse reaction
            tot_prod_rate += np.array(reverse_rate)
        elif diag_sp in products:
            rates['Production'][rea] = forward_rate # it gets produced in the forward reaction
            tot_prod_rate += np.array(forward_rate)
            rates['Destruction'][rev_rea] = reverse_rate # it is destroyed when reacting with something
            tot_dest_rate += np.array(reverse_rate)
    rates['Production']['total'] = tot_prod_rate
    rates['Production']['xlim_min'] = np.min(tot_prod_rate)*1e-1
    rates['Production']['xlim_max'] = np.max(tot_prod_rate)*5.
    rates['Destruction']['total'] = tot_dest_rate
    rates['Destruction']['xlim_min'] = np.min(tot_dest_rate)*1e-1
    rates['Destruction']['xlim_max'] = np.max(tot_dest_rate)*5.
    for k,v in rates['Production'].items():
        if any(v/tot_prod_rate > 1e-1): # use only >10% reactions
            important_rates['Production'].append(k)
    for k,v in rates['Destruction'].items():
        if any(v/tot_dest_rate > 1e-1):
            important_rates['Destruction'].append(k)
    return rates, important_rates
    
def get_prod_dest_reactions_to_plot(list_of_dat_list, diag_sp):
    reactions_to_plot = {'Production': ['total', 'xlim_min', 'xlim_max'], 'Destruction': ['total', 'xlim_min', 'xlim_max']}
    xlim_min = {'Production': 1e-1, 'Destruction': 1e-1}
    xlim_max = {'Production': 1e3, 'Destruction': 1e3}
    net_min, net_max = np.inf, -np.inf # to be able to minimise/maximise
    for dat_list in list_of_dat_list:
        for d in dat_list:
            r, imp = get_prod_dest_rates(d, diag_sp)
            reactions_to_plot['Production'] += imp['Production']
            reactions_to_plot['Destruction'] += imp['Destruction']
            xlim_min['Production'] = min(xlim_min['Production'], r['Production']['xlim_min'])
            xlim_min['Destruction'] = min(xlim_min['Destruction'], r['Destruction']['xlim_min'])
            xlim_max['Production'] = max(xlim_max['Production'], r['Production']['xlim_max'])
            xlim_max['Destruction'] = max(xlim_max['Destruction'], r['Destruction']['xlim_max'])
            net_min = min(net_min, np.min(r['Production']['total']-r['Destruction']['total']))
            net_max = max(net_max, np.max(r['Production']['total']-r['Destruction']['total']))
    reactions_to_plot['Production'] = list(set(reactions_to_plot['Production']))
    reactions_to_plot['Destruction'] = list(set(reactions_to_plot['Destruction']))
    return reactions_to_plot, xlim_min, xlim_max, net_min, net_max

def plot_prod_dest_rates(dat_list, param_list, diag_sp, rplot, xlim_lower, xlim_upper, figsave, sim_type):
    for i,d in enumerate(dat_list):
        prod_dest, _ = get_prod_dest_rates(d, diag_sp)
        prod_dest['Production'] = {key: prod_dest['Production'][key] for key in prod_dest['Production'] if key in rplot['Production']}
        prod_dest['Destruction'] = {key: prod_dest['Destruction'][key] for key in prod_dest['Destruction'] if key in rplot['Destruction']}
        labels = [r'k$_{prod}$ [cm$^{-3}$s$^{-1}$]', r'k$_{dest}$ [cm$^{-3}$s$^{-1}$]']
        fig, ax = plt.subplots(ncols=1, nrows=2, sharex = True, tight_layout = True)
        ax = ax.flatten()
        for j,rate_type in enumerate(['Production', 'Destruction']):
            for k,v in prod_dest[rate_type].items():
                if k in ['xlim_min', 'xlim_max']:
                    continue
                if k == 'total':
                    ax[j].plot(v, d['atm']['pco']/1e6, label = 'Total', c = 'k')
                else:
                    ax[j].plot(v, d['atm']['pco']/1e6, label = k)
            ax[j].set_yscale('log')
            ax[j].set_xscale('log')
            ax[j].set_ylabel('Pressure [bar]')
            ax[j].set_xlabel(labels[j])
            ax[j].set_xlim(xlim_lower[rate_type],xlim_upper[rate_type]) # more clever way to set this? by 1% of total rate or something like that?
            ax[j].invert_yaxis()
            #ax[j].legend(loc = 'upper left')
            ax[j].legend(bbox_to_anchor = (1,0.85))
        ax[0].set_title(legend_lab[sim_type].format(param_list[i]))
        if figsave:
            fig.savefig(plot_folder + 'prod_dest/detailed'+end_str[sim_type]+'_'+sim_names[i]+network+'.pdf', bbox_inches = 'tight')

def plot_prod_dest_rates_normed(dat_list, param_list, diag_sp, rplot, net_lower, net_upper, figsave, sim_type):
    for i,d in enumerate(dat_list):
        prod_dest, _ = get_prod_dest_rates(d, diag_sp)
        prod_dest['Production'] = {key: prod_dest['Production'][key] for key in prod_dest['Production'] if key in rplot['Production']}
        prod_dest['Destruction'] = {key: prod_dest['Destruction'][key] for key in prod_dest['Destruction'] if key in rplot['Destruction']}
        labels = [r'$X_{prod}$', r'$X_{dest}$', r'$k_{net}$ [cm$^{-3}$s$^{-1}$]']
        fig, ax = plt.subplots(ncols=1, nrows=3, sharey = True, tight_layout = True)
        ax = ax.flatten()
        for j,rate_type in enumerate(['Production', 'Destruction', 'total']):
            if rate_type != 'total':
                for k,v in prod_dest[rate_type].items():
                    if k in ['total', 'xlim_min', 'xlim_max']:
                        continue
                    else:
                        ax[j].plot(v/prod_dest[rate_type]['total'], d['atm']['pco']/1e6, label = k)
                ax[j].legend(bbox_to_anchor = (1,0.85))
                ax[j].set_xlim(-0.02,1.02)
            else:
                ax[j].plot(prod_dest['Production']['total']-prod_dest['Destruction']['total'], d['atm']['pco']/1e6, label = 'Total', c = 'k')
                ax[j].set_xlim(net_lower, net_upper)
                ax[j].set_xscale('symlog')
                ax[j].axvline(0, color = 'r', linestyle = '--')
            ax[j].set_yscale('log')
            ax[j].set_ylabel('Pressure [bar]')
            ax[j].set_xlabel(labels[j])
            ax[j].invert_yaxis()
        ax[0].set_title(legend_lab[sim_type].format(param_list[i]))
        if figsave:
            fig.savefig(plot_folder + 'prod_dest/normed_detailed'+end_str[sim_type]+'_'+sim_names[i]+network+'.pdf', bbox_inches = 'tight')
            
def plot_prod_dest_rates_normed_selected(dat_list, param_list, diag_sp, rplot, figsave, sim_type):
    legend_xanchors = [0.5, 0.5] # otherwise legends are all over the place, not sure why...
    legend_yanchors = [0.666, 0.323]
    prod_dest_archean, _ = get_prod_dest_rates(data_archean, diag_sp)
    fig, ax = plt.subplots(ncols=len(plot_idx[sim_type]), nrows=3, sharey = True, tight_layout = True)
    for i,idx in enumerate(plot_idx[sim_type]):
        if idx == 'archean': # if it is the archean simulation, use the archean data
            prod_dest = {'Production': {}, 'Destruction': {}}
            prod_dest['Production'] = {key: prod_dest_archean['Production'][key] for key in prod_dest_archean['Production'] if key in rplot['Production']}
            prod_dest['Destruction'] = {key: prod_dest_archean['Destruction'][key] for key in prod_dest_archean['Destruction'] if key in rplot['Destruction']}
            d = data_archean
        else:
            prod_dest, _ = get_prod_dest_rates(dat_list[idx], diag_sp)
            prod_dest['Production'] = {key: prod_dest['Production'][key] for key in prod_dest['Production'] if key in rplot['Production']}
            prod_dest['Destruction'] = {key: prod_dest['Destruction'][key] for key in prod_dest['Destruction'] if key in rplot['Destruction']}
            d = dat_list[idx]
        labels = [r'$X_{prod}$', r'$X_{dest}$', r'$k_{tot}$ [cm$^{-3}$s$^{-1}$]']
        for j,rate_type in enumerate(['Production', 'Destruction', 'total']):
            if rate_type != 'total':
                for k,v in prod_dest[rate_type].items():
                    if k in ['total', 'xlim_min', 'xlim_max']:
                        continue
                    else:
                        ax[j,i].plot(v/prod_dest[rate_type]['total'], d['atm']['pco']/1e6, label = k)
                if j in [0,1] and i == 0:
                    handles, leg_labels = ax[j,i].get_legend_handles_labels()
                    fig.tight_layout(h_pad = 4.2) # use if legends are below subplots
                    fig.legend(handles, leg_labels, loc = 'center', bbox_to_anchor = (legend_xanchors[j], legend_yanchors[j]), ncol = 3)
                ax[j,i].set_xlim(-0.02,1.02)
            else:
                ax[j,i].plot(prod_dest['Production']['total'], d['atm']['pco']/1e6, label = 'Production', c = 'r', ls = '-')
                ax[j,i].plot(prod_dest['Destruction']['total'], d['atm']['pco']/1e6, label = 'Destruction', c = 'k', ls = '-')
                ax[j,i].set_xlim(3e-5, 1e4)
                ax[j,i].set_xscale('log')
                ax[j,i].axvline(0, color = 'r', linestyle = '--')
                if i == 0:
                    ax[j,i].legend(loc = 'lower right')
            ax[j,i].set_yscale('log')
            ax[j,0].set_ylabel('Pressure [bar]')
            ax[j,i].set_xlabel(labels[j])
            ax[j,i].text(0.03, 0.93, '{})'.format(chr(97+j*4+i%4)), transform=ax[j,i].transAxes)
        if idx == 'archean': # add (Archean) to the title if it is the Archean simulation
            ax[0,i].set_title(legend_lab[sim_type].format(archean_params[sim_type]) + ' (Archean)')       
        else:
            ax[0,i].set_title(legend_lab[sim_type].format(param_list[idx]))
    ax[0,0].invert_yaxis() # due to sharey, only need to reverse one
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/selected_normed_detailed'+end_str[sim_type]+network+'.pdf', bbox_inches = 'tight')

def plot_prod_dest_rates_archean(dat, diag_sp, rplot, xlim_lower, xlim_upper, figsave):
    prod_dest, _ = get_prod_dest_rates(dat, diag_sp)
    prod_dest['Production'] = {key: prod_dest['Production'][key] for key in prod_dest['Production'] if key in rplot['Production']+['xlim_min', 'xlim_max']}
    prod_dest['Destruction'] = {key: prod_dest['Destruction'][key] for key in prod_dest['Destruction'] if key in rplot['Destruction']+['xlim_min', 'xlim_max']}
    labels = [r'k$_{prod}$ [cm$^{-3}$s$^{-1}$]', r'k$_{dest}$ [cm$^{-3}$s$^{-1}$]']
    fig, ax = plt.subplots(ncols=1, nrows=2, sharex = True, tight_layout = True)
    ax = ax.flatten()
    for j,rate_type in enumerate(['Production', 'Destruction']):
        for k,v in prod_dest[rate_type].items():
            if k in ['xlim_min', 'xlim_max']:
                continue
            if k == 'total':
                ax[j].plot(v, dat['atm']['pco']/1e6, label = 'Total', c = 'k')
            else:
                ax[j].plot(v, dat['atm']['pco']/1e6, label = k)
        ax[j].set_yscale('log')
        ax[j].set_xscale('log')
        ax[j].set_ylabel('Pressure [bar]')
        ax[j].set_xlabel(labels[j])
        ax[j].set_xlim(xlim_lower[rate_type],xlim_upper[rate_type])
        ax[j].invert_yaxis()
        ax[j].legend(bbox_to_anchor = (1,0.85))
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/archean'+network+'.pdf', bbox_inches = 'tight')

def plot_prod_dest_rates_archean_normed(dat, diag_sp, rplot, net_lower, net_upper, figsave):
    prod_dest, _ = get_prod_dest_rates(dat, diag_sp)
    prod_dest['Production'] = {key: prod_dest['Production'][key] for key in prod_dest['Production'] if key in rplot['Production']+['xlim_min', 'xlim_max']}
    prod_dest['Destruction'] = {key: prod_dest['Destruction'][key] for key in prod_dest['Destruction'] if key in rplot['Destruction']+['xlim_min', 'xlim_max']}
    labels = [r'$X_{prod}$', r'$X_{dest}$', r'$k_{net}$ [cm$^{-3}$s$^{-1}$]']
    fig, ax = plt.subplots(ncols=1, nrows=3, sharey = True, tight_layout = True)
    ax = ax.flatten()
    for j,rate_type in enumerate(['Production', 'Destruction', 'total']):
        if rate_type != 'total':
            for k,v in prod_dest[rate_type].items():
                if k in ['total', 'xlim_min', 'xlim_max']:
                    continue
                else:
                    ax[j].plot(v/prod_dest[rate_type]['total'], dat['atm']['pco']/1e6, label = k)
            ax[j].legend(bbox_to_anchor = (1,0.85))
            ax[j].set_xlim(-0.02,1.02)
        else:
            ax[j].plot(prod_dest['Production']['total']-prod_dest['Destruction']['total'], dat['atm']['pco']/1e6, label = 'Total', c = 'k')
            ax[j].set_xlim(net_lower, net_upper)
            ax[j].set_xscale('symlog')
            ax[j].axvline(0, color = 'r', linestyle = '--')
        ax[j].set_yscale('log')
        ax[j].set_ylabel('Pressure [bar]')
        ax[j].set_xlabel(labels[j])
        ax[j].invert_yaxis()
    if figsave:
        fig.savefig(plot_folder + 'prod_dest/normed_archean'+network+'.pdf', bbox_inches = 'tight')

def plot_pt(dat_list, param_list, sim_type, figsave):
    fig, ax = plt.subplots(tight_layout = True)
    for d,p in zip(dat_list, param_list):
        ax.plot(d['atm']['Tco'], d['atm']['pco']/1e6, label = legend_lab[sim_type].format(p))
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylabel('Pressure [bar]')
    ax.set_xlabel('T [K]')  
    ax.legend(bbox_to_anchor = (1,0.95))
    if figsave:
        fig.savefig(plot_folder + 'TPs/'+end_str[sim_type][1:]+'.pdf', bbox_inches = 'tight')

def plot_rainrates_hcn_watercon_air_PT(list_of_dat_lists, list_of_param_lists, list_of_hcn_rain_lists, figsave):
    fig, ax = plt.subplots(nrows = 4, ncols = 4, figsize = (24,27), tight_layout = True)#, figsize = (22,18)) # take out tight layout here if legends are below subplots
    ax = ax.flatten()
    sim_types = ['BC', 'CtoO', 'dist', 'star']
    # legends below subplots
    legend_xanchors = [0.5, 0.5, 0.5, 0.5] # otherwise legends are all over the place, not sure why...
    legend_yanchors = [0.756, 0.503, 0.243, -0.015]
    hcn_rain_archean = rainout(data_archean)
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    i = 0
    for dat_list,param_list,hcn_rain_list in zip(list_of_dat_lists, list_of_param_lists, list_of_hcn_rain_lists):
        st = sim_types[i//4]
        # plotting hcn rain rates in zeroth column
        for p,r,c in zip(param_list, hcn_rain_list, colours):
            ax[i+0].plot(p, r, color = c, linestyle = '', marker = 'o', markersize = 10)
        ax[i+0].plot(archean_params[st], hcn_rain_archean, color = archean_colour, marker = archean_marker, markersize = 10)
        ax[i+0].set_ylabel(r'HCN rain-out rate [kg m$^{-2}$ yr$^{-1}$]')
        if i != 0:
            ax[i+0].set_yscale('log')
        ax[i+0].set_xscale(xscale[st])
        ax[i+0].set_xlabel(xlab[st])
        ax[i+0].text(0.03, 0.93, '{})'.format(chr(97+i+0)), transform=ax[i+0].transAxes)
        # plotting HCN and condensed water vertical structure and P-T profile (only star and dist simulations) in first, second and third columns, respectively
        for d,p in zip(dat_list,param_list):
            ax[i+1].plot(d['variable']['ymix'][:, d['variable']['species'].index('HCN')], d['atm']['pco']/1e6, label = legend_lab[st].format(p))
            ax[i+2].plot(d['variable']['ymix'][:, d['variable']['species'].index('H2O_l_s')], d['atm']['pco']/1e6)
            if sim_types[i//4] == 'star' or sim_types[i//4] == 'dist':
                ax[i+3].plot(d['atm']['Tco'], d['atm']['pco']/1e6)
        ax[i+1].set_xscale('log')
        ax[i+1].set_xlabel('X(HCN)')
        ax[i+1].set_yscale('log')
        ax[i+1].set_ylabel('Pressure [bar]')
        ax[i+1].invert_yaxis()
        ax[i+1].set_xlim((1e-15,1e-2))
        ax[i+1].text(0.03, 0.93, '{})'.format(chr(97+i+1)), transform=ax[i+1].transAxes)
        #ax[i+1].legend()
        ax[i+2].set_xscale('log')
        ax[i+2].set_xlabel('X(cloud)')
        ax[i+2].set_xlim((1e-15,1e-2)) # show relevant vertical mixing ratio interval
        ax[i+2].set_yscale('log')
        ax[i+2].set_ylabel('Pressure [bar]')
        ax[i+2].set_ylim((1e-4,1e0)) # zoom on relevant pressure interval for clouds
        ax[i+2].invert_yaxis()
        ax[i+2].text(0.03, 0.93, '{})'.format(chr(97+i+2)), transform=ax[i+2].transAxes)
        # plotting the fixed T-P profiles for BC and CtoO simulations and setting scales and labels for all
        if sim_types[i//4] == 'BC' or sim_types[i//4] == 'CtoO':
            ax[i+3].plot(data_archean['atm']['Tco'], data_archean['atm']['pco']/1e6)
        ax[i+3].set_xlabel('T [K]')
        ax[i+3].set_yscale('log')
        ax[i+3].set_ylabel('Pressure [bar]')
        ax[i+3].invert_yaxis()
        ax[i+3].set_xlim((125,385))
        ax[i+3].text(0.03, 0.93, '{})'.format(chr(97+i+3)), transform=ax[i+3].transAxes)
        handles, labels = ax[i+1].get_legend_handles_labels()
        fig.tight_layout(h_pad = 5.8) # use if legends are below subplots
        fig.legend(handles, labels, loc = 'center', bbox_to_anchor = (legend_xanchors[i//4], legend_yanchors[i//4]), ncol = 5)#, ncols = 2) took out for legends on side, also use upper right for loc
        i += 4
    
    if figsave:
        fig.savefig(plot_folder + 'rainout_rates/rain_vertical_pt'+network+'.pdf', bbox_inches = 'tight')
    
def get_rad_prof(star):
    ''' Taken from parallel_functions.py but changed so relative passes are correct. 
        Gives back the location of the needed stellar radiation profile file for a given star.'''
    rad_file = ''
    if star == 'SUN':
        rad_file = '../atm/stellar_flux/Gueymard_solar.txt'
    else:
        rad_file = '/scratch/s2555875/stellar_flux/' + star.lower() + '.txt' # just to make sure it is lower case
    return rad_file

def plot_stellar_spectra(figsave):
    fig, ax = plt.subplots(nrows = 3, ncols = 1, sharex = True, sharey = True, tight_layout = True, figsize = (8, 10))
    ax = ax.flatten()
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i,star_name in enumerate(list(star_df.Name)):
        if star_name == 'SUN':
            star = np.genfromtxt('/home/s2555875/VULCAN-2/atm/stellar_flux/Gueymard_solar.txt', comments = '#', names = ['lambda', 'flux'])
        else:
            star = np.genfromtxt(stellar_spectra_folder + '/{}.txt'.format(star_name).lower(), comments = '#', names = ['lambda', 'flux'])
        ax[i//5].plot(star['lambda'], star['flux'], label = star_name, c = colours[i])
    ax[2].set_xlabel('Wavelength [nm]')
    ax[1].set_ylabel(r'Flux [ergs/cm$^2$/s/nm]')
    for i in range(3):
        ax[i].set_yscale('log')
        ax[i].set_xscale('log')
    ax[0].set_xlim((2, 700)) # VULCAN uses this range of the spectrum for photochemistry
    ax[0].set_ylim((9e0, 3e8))
    fig.legend(bbox_to_anchor=(1.22,0.72))
    
    if figsave:
        fig.savefig(plot_folder + 'spectra_comp/stellar_spectra_comp.pdf', bbox_inches = 'tight')
#%%
# boundary condition case
data_bc, bc_flux = read_in('BC', nsim)

hcn_rain, rain = [], [] # storing the rainout rates of HCN and water
for d in data_bc:
    hcn_rain.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    rain.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

plot_vertical_n(data_bc, 'HCN', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'H2O_l_s', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'HNCO', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'H2CN', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'C2H3CN', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'C2H3', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'C2H6', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'CH4', bomb_rate, 'BC', figsave = True)
plot_vertical_n(data_bc, 'CH3', bomb_rate, 'BC', figsave = True)
plot_end_time(data_bc, figsave = True, sim_type = 'BC')
plot_evo_layer(data_bc, bomb_rate, 'HCN', 0, figsave = True, sim_type = 'BC')
plot_convergence(data_bc, figsave = True, sim_type = 'BC')
plot_rain(hcn_rain, bomb_rate, 'BC', extra_list = bc_flux, rain_spec = 'HCN_rain', figsave = True)
plot_rain(rain, bomb_rate, 'BC', extra_list = bc_flux, rain_spec = 'H2O_rain', figsave = True)
plot_tot_rate(data_bc, bomb_rate, 'BC', diag_sp = 'HCN', figsave = True)
plot_prod_dest(data_bc, bomb_rate, 'BC', diag_sp = 'HCN', figsave = True)
plot_prod_dest_layer(data_bc, bomb_rate, 'BC', 0, diag_sp = 'HCN', figsave = True)
plot_prod_dest_many_layer(data_bc, bomb_rate, 'BC', pressure_levels, 4, diag_sp = 'HCN', figsave = True)
#%%
# C/O case
data_CtoO, C_to_O = read_in('CtoO', nsim)

hcn_rain_CtoO, rain_CtoO = [], []
for d in data_CtoO:
    hcn_rain_CtoO.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    rain_CtoO.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
plot_vertical_n(data_CtoO, 'HCN', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'H2O_l_s', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'HNCO', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'H2CN', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'C2H3CN', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'C2H3', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'C2H6', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'CH4', C_to_O, 'CtoO', figsave = True)
plot_vertical_n(data_CtoO, 'CH3', C_to_O, 'CtoO', figsave = True)
plot_vertical_many(data_CtoO, C_to_O, 'CtoO', True)
plot_end_time(data_CtoO, figsave = True, sim_type = 'CtoO')
#plot_evo_layer(data_CtoO, C_to_O, 'HCN', 0, figsave = True)
plot_convergence(data_CtoO, figsave = True, sim_type = 'CtoO')
plot_rain(hcn_rain_CtoO, C_to_O, 'CtoO', figsave = True, rain_spec = 'HCN_rain')
plot_rain(rain_CtoO, C_to_O, 'CtoO', figsave = True, rain_spec = 'H2O_rain')
plot_tot_rate(data_CtoO, C_to_O, 'CtoO', diag_sp = 'HCN', figsave = True)
plot_tot_rate_selected(data_CtoO, C_to_O, 'CtoO', diag_sp = 'HCN', figsave = True)
plot_prod_dest(data_CtoO, C_to_O, 'CtoO', diag_sp = 'HCN', figsave = True)
plot_prod_dest_layer(data_CtoO, C_to_O, 'CtoO', 0, diag_sp = 'HCN', figsave = True)
plot_prod_dest_many_layer(data_CtoO, C_to_O, 'CtoO', pressure_levels, 4, diag_sp = 'HCN', figsave = True)
plot_prod_dest_selected(data_CtoO, C_to_O, 'CtoO', 'HCN', True)
# %%
# star case
data_star = read_in('star', number_of_sim = 13)
hcn_rain_star, rain_star = [], []

for d in data_star:
    hcn_rain_star.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    rain_star.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
plot_vertical_n(data_star, 'HCN', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'H2O_l_s', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'HNCO', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'H2CN', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'C2H3CN', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'C2H3', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'C2H6', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'CH4', T_eff, 'star', figsave = True)
plot_vertical_n(data_star, 'CH3', T_eff, 'star', figsave = True)
plot_vertical_many(data_star, T_eff, 'star', True)
plot_end_time(data_star, figsave = True, sim_type = 'star')
plot_evo_layer(data_star, T_eff, 'HCN', 0, figsave = True, sim_type = 'star')
plot_convergence(data_star, figsave = True, sim_type = 'star')
plot_rain(hcn_rain_star, T_eff, 'star', figsave = True, rain_spec = 'HCN_rain')
plot_rain(rain_star, T_eff, 'star', figsave = True, rain_spec = 'H2O_rain')
plot_pt(data_star, T_eff, 'star', figsave = True)
plot_tot_rate(data_star, T_eff, 'star', diag_sp = 'HCN', figsave = True)
plot_tot_rate_selected(data_star, T_eff, 'star', diag_sp = 'HCN', figsave = True)
plot_prod_dest(data_star, T_eff, 'star', diag_sp = 'HCN', figsave = True)
plot_prod_dest_layer(data_star, T_eff, 'star', 0, diag_sp = 'HCN', figsave = True)
plot_prod_dest_many_layer(data_star, T_eff, 'star', pressure_levels, 4, diag_sp = 'HCN', figsave = True)
plot_prod_dest_selected(data_star, T_eff, 'star', 'HCN', True)
plot_stellar_spectra(figsave = True)
# %%
# distance case
data_dist, T_surf = read_in('dist', nsim)
hcn_rain_dist, rain_dist = [], []

for d in data_dist:
    hcn_rain_dist.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    rain_dist.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
plot_vertical_n(data_dist, 'HCN', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'H2O_l_s', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'HNCO', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'H2CN', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'C2H3CN', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'C2H3', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'C2H6', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'CH4', a_list, 'dist', figsave = True)
plot_vertical_n(data_dist, 'CH3', a_list, 'dist', figsave = True)
plot_vertical_many(data_dist, a_list, 'dist', True)
plot_end_time(data_dist, figsave = True, sim_type = 'dist')
plot_evo_layer(data_dist, a_list, 'HCN', 0, figsave = True, sim_type = 'dist')
plot_convergence(data_dist, figsave = True, sim_type = 'dist')
plot_rain(hcn_rain_dist, a_list, 'dist', extra_list = T_surf, figsave = True, rain_spec = 'HCN_rain')
plot_rain(rain_dist, a_list, 'dist', extra_list = T_surf, figsave = True, rain_spec = 'H2O_rain')
plot_pt(data_dist, a_list, 'dist', figsave = True)
plot_tot_rate(data_dist, a_list, 'dist', diag_sp = 'HCN', figsave = True)
plot_tot_rate_selected(data_dist, a_list, 'dist', diag_sp = 'HCN', figsave = True)
plot_prod_dest(data_dist, a_list, 'dist', diag_sp = 'HCN', figsave = True)
plot_prod_dest_layer(data_dist, a_list, 'dist', 0, diag_sp = 'HCN', figsave = True)
plot_prod_dest_many_layer(data_dist, a_list, 'dist', pressure_levels, 4, diag_sp = 'HCN', figsave = True)
plot_prod_dest_selected(data_dist, a_list, 'dist', 'HCN', True)
#%%
# reaction rate plots
pr.reset_plt(ticksize = 15, fontsize = 17, fxsize = 11, fysize = 10)
r_to_lot, xmin, xmax, netmin, netmax = get_prod_dest_reactions_to_plot([data_bc, data_CtoO, data_dist, data_star], 'HCN')
# boundary condition case
plot_prod_dest_rates(data_bc, bomb_rate, 'HCN', r_to_lot, xmin, xmax, figsave = True, sim_type = 'BC')
# C/O case
plot_prod_dest_rates(data_CtoO, C_to_O, 'HCN', r_to_lot, xmin, xmax, figsave = True, sim_type = 'CtoO')
# distance case
plot_prod_dest_rates(data_dist, a_list, 'HCN', r_to_lot, xmin, xmax, figsave = True, sim_type = 'dist')
# star case
plot_prod_dest_rates(data_star, T_eff, 'HCN', r_to_lot, xmin, xmax, figsave = True, sim_type = 'star')
# archean case
plot_prod_dest_rates_archean(data_archean, 'HCN', r_to_lot, xmin, xmax, figsave = True)
#%%
pr.reset_plt(ticksize = 15, fontsize = 17, fxsize = 11, fysize = 15)
plot_prod_dest_rates_normed(data_bc, bomb_rate, 'HCN', r_to_lot, netmin, netmax, figsave = True, sim_type = 'BC')
plot_prod_dest_rates_normed(data_CtoO, C_to_O, 'HCN', r_to_lot, netmin, netmax, figsave = True, sim_type = 'CtoO')
plot_prod_dest_rates_normed(data_dist, a_list, 'HCN', r_to_lot, netmin, netmax, figsave = True, sim_type = 'dist')
plot_prod_dest_rates_normed(data_star, T_eff, 'HCN', r_to_lot, netmin, netmax, figsave = True, sim_type = 'star')
plot_prod_dest_rates_archean_normed(data_archean, 'HCN', r_to_lot, netmin, netmax, figsave = True)
#%%
pr.reset_plt(ticksize = 15, fontsize = 17, fxsize = 24, fysize = 18)
plot_prod_dest_rates_normed_selected(data_CtoO, C_to_O, 'HCN', r_to_lot, figsave = True, sim_type = 'CtoO')
plot_prod_dest_rates_normed_selected(data_dist, a_list, 'HCN', r_to_lot, figsave = True, sim_type = 'dist')
plot_prod_dest_rates_normed_selected(data_star, T_eff, 'HCN', r_to_lot, figsave = True, sim_type = 'star')
#%%
pr.reset_plt(ticksize = 16, fontsize = 19, fxsize = 24, fysize = 27)
plot_rainrates_hcn_watercon_air_PT([data_bc, data_CtoO, data_dist, data_star], [bomb_rate, C_to_O, a_list, T_eff], [hcn_rain, hcn_rain_CtoO, hcn_rain_dist, hcn_rain_star], figsave = True)
#%%
