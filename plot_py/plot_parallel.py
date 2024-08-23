#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

# general parameters
nsim = 15 # hardcoded for now, change later...
out_folder = '/scratch/s2555875/output/'
plot_folder = '/scratch/s2555875/plot/'

# setting up the boundary condition case
bomb_rate = np.linspace(3e23, 1e25, nsim) # values from Pearce et al. (2022) Fig A4 min and max
bc_folder= '/scratch/s2555875/BC_files/'
bc_spec = 'H2'
bc_linestyle = '-'

# setting up the C/O case
co2_for_CtoO_range = np.linspace(0,0.1,nsim, endpoint = True)
mixing_folder = '/scratch/s2555875/mixing_files/'

# setting up the distance case
helios_output_folder = '/scratch/s2555875/HELIOS/output/'

# base simulation of Archean
with open(out_folder+'archean.vul', 'rb') as handle:
    data_archean = pickle.load(handle)
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

def read_in(sim_type, number_of_sim, start_number = '0', start_str = ''):
    ''' Reads in all the data for a specific simulation type. The simulation type dictates how the files
        are found and what extra parameters this function returns. It builds on the idea that output
        files have a format of start_str+sim_+(start_number+)sim_number+sim_type+.vul.'''
    dat_list = [] # list to store the results (dicts)
    H2_flux = []
    T_surface = []
    ctoo = []
    extra_str = '' # to get the proper output
    if sim_type == 'BC':
        #extra_str = '_onlyH2.vul'
        extra_str = '_meteor.vul'
    elif sim_type == 'CtoO':
        extra_str = '_CtoO.vul'
    elif sim_type == 'star':
        extra_str = '_star.vul'
    elif sim_type == 'dist':
        extra_str = '_dist.vul'
    for i in range(number_of_sim):
        sim = start_str + 'sim_' # setting up the names of files
        if i < 10:
            sim += start_number + str(i)
        else:
            sim += str(i)
        sim += extra_str

        with open(out_folder + sim, 'rb') as handle:
            data_i = pickle.load(handle)
            dat_list.append(data_i)

        if sim_type == 'BC':
            with open(bc_folder+'BC_bot_'+sim[:-4]+'.txt') as f: # vas sim[:-11]+'.txt
                for line in f:
                    lin = line.split()
                    if line[0] != '#' and lin[0] == bc_spec:
                        H2_flux.append(float(lin[1]))
                        break
        elif sim_type == 'CtoO':
            ctoo.append(calc_C_to_O(data_i, mixing_folder + sim[:-4] + 'mixing.txt'))
        elif sim_type == 'dist':
            surface_temp = np.genfromtxt(helios_output_folder + sim[:-4] + '/{}_tp.dat'.format(sim[:-4]), dtype = None, skip_header = 2, usecols = (1))[0]
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
    rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = rain_rate * 5.237e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)

def plot_rain(hcn_rain_list, param_list, sim_type, figname = None, f_size = 13, l_size = 12, bc_flux_list = [], surf_temp = [], yscale = 'log', mol = 'HCN', plot_Pearce = True):
    ''' Function to plot rainout rates. Specific cases for different simulation types are included via
        conditions. If needed, using the plot_Pearce parameter, it plots the fiducial model value which
        is marked with a B at the end of the variables.'''
    if plot_Pearce:
        with open(out_folder+'B_nofix.vul', 'rb') as handle: # for comparison
            data_B = pickle.load(handle)
        hcn_rain_B = 3.4e-12 #rainout(data_B)
        bc_B = 2.3e+10
        bomb_B = 1.2e24 * 2.3/3.42 # maybe it is not scaled or my calc is bad
        C_to_O_B = calc_C_to_O(data_B, '../atm/mixing_Pearce_B.txt')
    colour_b = 'tab:blue'
    fig, ax = plt.subplots(tight_layout = True)
    if plot_Pearce:
        ax.plot(param_list, hcn_rain_list, color = colour_b, label = 'own')
    else:
        ax.plot(param_list, hcn_rain_list, color = colour_b)
    if sim_type == 'BC': # plotting comparison and H2 flux
        if plot_Pearce:
            ax.plot(bomb_B, hcn_rain_B, marker = '*', linestyle = '', label = 'Pearce et al. (2022)')
        ax.set_xscale('log')
        ax.set_xlabel(r'Mass delivery rate [g Gyr$^{-1}$]', fontsize = f_size)
        ax.tick_params(axis = 'y', labelcolor = colour_b, labelsize = l_size)
        ax.set_ylabel(mol + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]', color = colour_b, fontsize = f_size)
        ax.tick_params(which='both', direction='out', width=1, length = 4)
        ax.tick_params(axis = 'both', labelsize = l_size)
        ax1 = ax.twinx()
        colour = 'tab:orange'
        ax1.plot(param_list, bc_flux_list, color = colour, linestyle = bc_linestyle)
        if plot_Pearce:
            ax1.plot(bomb_B, bc_B, marker = '*', color = colour, linestyle = bc_linestyle, label = bc_spec)
        ax1.set_ylabel(r'H2 flux from surface [cm$^{-2}$ s$^{-1}$]', color = colour, fontsize = f_size)
        ax1.set_yscale('log')
        #ax1.set_xscale('log')
        #ax1.legend(loc = 'center right')
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(l_size)
        ax1.tick_params(which = 'both', direction = 'out', width = 1, length = 4)
        ax1.tick_params(axis = 'both', labelsize = l_size)
        ax1.tick_params(axis = 'y', labelcolor = colour)

    elif sim_type == 'C_to_O': # plotting comparison
        if plot_Pearce:
            ax.plot(C_to_O_B, hcn_rain_B, linestyle = '', marker = '*', color = colour_b, label = 'Pearce et al. (2022)')
        ax.set_xlabel('C/O', fontsize = f_size)
        ax.set_ylabel(mol + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]', fontsize = f_size)
        ax.tick_params(which='both', direction='out', width=1, length = 4)
        ax.tick_params(axis = 'both', labelsize = l_size)

    elif sim_type == 'star':
        if plot_Pearce:
            ax.plot(5332, hcn_rain_B, linestyle = '', marker = '*', color = colour_b, label = 'Pearce et al. (2022)')
        ax.set_xlabel(r'T$_{eff}$ [K]', fontsize = f_size)
        ax.set_ylabel(mol + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]', fontsize = f_size)
        ax.tick_params(which='both', direction='out', width=1, length = 4)
        ax.tick_params(axis = 'both', labelsize = l_size)
    
    elif sim_type == 'dist':
        if plot_Pearce:
            ax.plot(1, hcn_rain_B, marker = '*', linestyle = '', label = 'Pearce et al. (2022)')
        ax.set_xlabel('Distance [AU]', fontsize = f_size)
        ax.set_ylabel(mol + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]', fontsize = f_size)
        ax.tick_params(which='both', direction='out', width=1, length = 4)
        ax.tick_params(axis = 'both', labelsize = l_size)
        ax1 = ax.twiny()
        ax1.plot(surf_temp, hcn_rain_list, linestyle = '')
        ax1.set_xlabel(r'T$_{surf}$ [K]')
        ax1.invert_xaxis()
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(l_size)
    
    ax.set_yscale(yscale)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(l_size) 
    ax.legend()

    if figname != None:
        fig.savefig(plot_folder + figname)

def calc_rate(dat, spec_list, re_id, n):
    ''' Calculates the reaction rate by multiplying the reaction coefficient with the number densty of
        the reagents.'''
    rate = dat['variable']['k'][re_id][n]
    for sp in spec_list:
        if sp != 'M':
            rate *= dat['variable']['y'][n, dat['variable']['species'].index(sp)]
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
    
def plot_vertical_n(dat_list, spec, param_list, sim_type, figname = None, f_size = 13, l_size = 12):
    ''' Plots the vertical profile of the given species (spec) at the end of simulation for a list
        of simulations.'''
    fig, ax = plt.subplots(tight_layout = True, figsize = (10,6))
    cm = plt.get_cmap('nipy_spectral')
    n_colours = len(dat_list)
    lstyles = ['-', '--', '-.', ':']
    ax.set_prop_cycle(color=[cm(1.*i/n_colours) for i in range(n_colours)])
    for i in range(len(dat_list)):
        if sim_type == 'BC':
            ax.plot(dat_list[i]['variable']['y'][:, dat_list[i]['variable']['species'].index(spec)], dat_list[i]['atm']['zco'][1:]/1e5, label = r'$\dot{M}_{del}$ = ' + '{:.2e} g/Gyr'.format(param_list[i]), linestyle = lstyles[i%len(lstyles)])
        elif sim_type == 'C_to_O':
            ax.plot(dat_list[i]['variable']['y'][:, dat_list[i]['variable']['species'].index(spec)], dat_list[i]['atm']['zco'][1:]/1e5, label = 'C/O = {:.2f}'.format(param_list[i]), linestyle = lstyles[i%len(lstyles)])
        elif sim_type == 'star':
            ax.plot(dat_list[i]['variable']['y'][:, dat_list[i]['variable']['species'].index(spec)], dat_list[i]['atm']['zco'][1:]/1e5, label = r'T$_{eff}$'+'= {:.2f} K'.format(param_list[i]), linestyle = lstyles[i%len(lstyles)])
        elif sim_type == 'dist':
            ax.plot(dat_list[i]['variable']['y'][:, dat_list[i]['variable']['species'].index(spec)], dat_list[i]['atm']['zco'][1:]/1e5, label = 'd = {:.2f} AU'.format(param_list[i]), linestyle = lstyles[i%len(lstyles)])
        
    ax.set_xscale('log')
    ax.set_xlabel(r'n [cm$^{-3}$]', fontsize = f_size)
    ax.set_ylabel('Height [km]', fontsize = f_size)
    ax.tick_params(which='both', direction='out', width=1, length = 4)
    ax.tick_params(axis = 'both', labelsize = l_size)
    ax.legend(bbox_to_anchor = (1.1, 0.95), fontsize = str(l_size-4))
    if figname != None:
        fig.savefig(plot_folder + figname, bbox_inches = 'tight')

def plot_end_time(dat_list, figname = None, f_size = 13, l_size = 12):
    ''' Plots the end-of-simulation times for a list of simulations. It is used as a way of 
        seeing what difference there are between both converged and not converged simulations.'''
    fig, ax = plt.subplots(tight_layout = True)
    for i in range(len(dat_list)):
        ax.plot(i, dat_list[i]['variable']['t'], linestyle = '', marker = 'o', color = 'red')
    ax.set_yscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Simulation number', fontsize = f_size)
    ax.set_ylabel('End of simulation time [s]', fontsize = f_size)
    ax.tick_params(which='both', direction='out', width=1, length = 4)
    ax.tick_params(axis = 'both', labelsize = l_size)
    if figname != None:
        fig.savefig(plot_folder + figname)

def plot_evo_layer(dat_list, spec, layer = 0, figname = None, f_size = 13, l_size = 12):
    ''' Plots the evolution of a given species in a given layer for a list of simulations.'''
    fig, ax = plt.subplots(figsize = (12, 8))
    for i in range(len(dat_list)):
        ax.plot(dat_list[i]['variable']['t_time'], dat_list[i]['variable']['y_time'][:, 0, dat_list[i]['variable']['species'].index(spec)], label = i)
    ax.set_xlabel('Time [s]', fontsize = f_size)
    ax.set_ylabel('n [cm-3]', fontsize = f_size)
    ax.set_yscale('log')
    ax.legend()
    ax.set_ylim(1e-2,None)
    ax.tick_params(which='both', direction='out', width=1, length = 4)
    ax.tick_params(axis = 'both', labelsize = l_size)
    if figname != None:
        fig.savefig(plot_folder + figname)

def plot_convergence(dat_list, figname = None, f_size = 13, l_size = 12):
    ''' It checks and plots whether the convergence criteria has been met in all simulations in the given list. 
        Template for calculation is taken from the VULCAN code (Tsai et al 2017, 2020).'''
    yconv_cri = 0.01
    yconv_min = 0.1
    slope_cri = 1.e-4
    fig, ax = plt.subplots(tight_layout = True)
    for i in range(len(dat_list)):
        longdy = dat_list[i]['variable']['longdy']
        longdydt = dat_list[i]['variable']['longdydt']
        slope_min = min( np.amin(dat_list[i]['atm']['Kzz']/(0.1*dat_list[i]['atm']['Hp'][:-1])**2) , 1.e-8)
        slope_min = max(slope_min, 1.e-10)
        if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min):
            ax.plot(i, 'Yes', 'ro')
        else:
            ax.plot(i, 'No', 'ro')
    ax.set_xlabel('Simulation number', fontsize = f_size)
    ax.set_ylabel('Converged', fontsize = f_size)
    ax.tick_params(which='both', direction='out', width=1, length = 4)
    ax.tick_params(axis = 'both', labelsize = l_size)
    
    if figname != None:
        fig.savefig(plot_folder + figname)

def plot_rain_converged(dat_list, rain_list, param_list, sim_type, figname = None, f_size = 15, l_size = 14, rain_spec = 'HCN_rain', extra_list = [], plot_non_conv = False):
    ''' Plots the rainout rates for a list of simulations for a given simulation types. It 
        distinguishes between converged and non-converged simulations (full and empty circles, respectively).
        Plotting and convergence caalculations are taken from previous functions.'''
    yconv_cri = 0.01
    yconv_min = 0.1
    slope_cri = 1.e-4
    conv_rain_list, conv_param_list, conv_extra_list = [], [], []
    non_conv_rain_list, non_conv_param_list, non_conv_extra_list = [], [], []
    for i in range(len(dat_list)):
        longdy = dat_list[i]['variable']['longdy']
        longdydt = dat_list[i]['variable']['longdydt']
        slope_min = min( np.amin(dat_list[i]['atm']['Kzz']/(0.1*dat_list[i]['atm']['Hp'][:-1])**2) , 1.e-8)
        slope_min = max(slope_min, 1.e-10)
        if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min):
            conv_rain_list.append(rain_list[i])
            conv_param_list.append(param_list[i])
            if extra_list:
                conv_extra_list.append(extra_list[i])
        else:
            non_conv_rain_list.append(rain_list[i])
            non_conv_param_list.append(param_list[i])
            if extra_list:
                non_conv_extra_list.append(extra_list[i])
    gpm = 0.
    if rain_spec == 'HCN_rain':
        gpm = 27
    elif rain_spec == 'H2O_rain':
        gpm = 18
    rain_archean = rainout(data_archean, rain_spec = rain_spec, g_per_mol = gpm)
    param_archean = 0.
    if sim_type == 'BC':
        param_archean = 1.2e24 * 2.3/3.42 # scaled bomb rate
    elif sim_type == 'CtoO':
        param_archean = calc_C_to_O(data_archean, '/home/s2555875/VULCAN-2/atm/mixing_table_archean.txt')
    elif sim_type == 'star':
        param_archean = 5600.
    elif sim_type == 'dist':
        param_archean = 1.

    popt, pcov = curve_fit(lin, conv_param_list, conv_rain_list)
    param_x = np.linspace(conv_param_list[0], conv_param_list[-1], 42, endpoint = True)
    fig, ax = plt.subplots(tight_layout = True)
    ax.plot(conv_param_list, conv_rain_list, color = 'navy', linestyle = '', marker = '.', markersize = 10)
    if plot_non_conv and non_conv_param_list: # told to plot non converged and there are non converged sims
        ax.plot(non_conv_param_list, non_conv_rain_list, color = 'navy', linestyle = '', marker = '.', fillstyle = 'none', markersize = 10)
    ax.plot(param_archean, rain_archean, color = 'darkorange', marker = '*', markersize = 10)
    ax.plot(param_x, lin(param_x, popt[0], popt[1]), linestyle = '--', alpha = 0.4, color = 'r')
    ax.set_yscale('log')
    if sim_type == 'BC':
        ax.set_xscale('log')
        ax.set_xlabel(r'Mass delivery rate [g Gyr$^{-1}$]', fontsize = f_size)
        ax1 = ax.twiny()
        ax1.plot(conv_extra_list, conv_rain_list, linestyle = '')
        if plot_non_conv and non_conv_param_list: # told to plot non converged and there are non converged sims
            ax1.plot(non_conv_extra_list, non_conv_rain_list, color = 'navy', linestyle = '', marker = '.', markerfacecolor = 'none')
        ax1.set_xlabel(r'H$_2$ flux [cm$^{-2}$ s$^{-1}$]', fontsize = f_size)
        ax1.set_xscale('log')
        ax1.tick_params(which = 'both', direction='out', width=1, length = 4)
        ax1.tick_params(axis = 'both', labelsize = l_size)
        #for tick in ax1.xaxis.get_major_ticks():
        #    tick.label1.set_fontsize(l_size)
    elif sim_type == 'CtoO':
        ax.set_xlabel('C/O', fontsize = f_size)
    elif sim_type == 'star':
        ax.set_xlabel(r'T$_{eff}$ [K]', fontsize = f_size)
    elif sim_type == 'dist':
        ax.set_xlabel('Distance [AU]', fontsize = f_size)
        ax1 = ax.twiny()
        ax1.plot(conv_extra_list, conv_rain_list, linestyle = '')
        if plot_non_conv and non_conv_param_list: # told to plot non converged and there are non converged sims
            ax1.plot(non_conv_extra_list, non_conv_rain_list, color = 'navy', linestyle = '', marker = '.', markerfacecolor = 'none')
        ax1.set_xlabel(r'T$_{surf}$ [K]', fontsize = f_size)
        ax1.invert_xaxis()
        ax1.tick_params(which = 'both', direction='out', width=1, length = 4)
        ax1.tick_params(axis = 'both', labelsize = l_size)
        

    ax.set_ylabel(rain_spec[:-5] + r' rain-out rate [kg m$^{-2}$ yr$^{-1}$]', fontsize = f_size)
    ax.tick_params(which='both', direction='out', width=1, length = 4)
    ax.tick_params(axis = 'both', labelsize = l_size)

    if figname != None:
        fig.savefig(plot_folder + figname)
    

#%%
# BC case
data_bc, bc_flux = read_in('BC', nsim)
#nsim_extra = 8
#data_bc_extra, bc_flux_extra = read_in('BC', nsim_extra, start_number = '00') # post simulations
#data_bc = [data_bc[0]] + data_bc_extra + data_bc[1:] # to have bombardment rates in order
#bc_flux = [bc_flux[0]] + bc_flux_extra + bc_flux[1:]
#bomb_rate_extra = np.linspace(3.5e23, 9e23, nsim_extra)
#bomb_rate = np.concatenate((np.concatenate((np.array([bomb_rate[0]]), bomb_rate_extra)), bomb_rate[1:]))

hcn_rain = [] # storing the rainout rates
for d in data_bc:
    hcn_rain.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))

rain = [] # storing the rainout rates
for d in data_bc:
    rain.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

plot_rain(rain, bomb_rate, 'BC', figname = 'H2O_rainout_meteor_onlyH2.pdf', bc_flux_list = bc_flux, mol = 'Water')
plot_rain(hcn_rain, bomb_rate, 'BC', figname = 'HCN_rainout_meteor_onlyH2.pdf', bc_flux_list = bc_flux)
plot_vertical_n(data_bc, 'HCN', bomb_rate, 'BC', figname = 'HCN_air_meteoritic.pdf')
plot_end_time(data_bc, figname = 'end_time_meteor.pdf')
plot_evo_layer(data_bc, 'HCN', figname = 'hcn_evo_meteor.pdf')
plot_convergence(data_bc, figname = 'convergence_meteor.pdf')
#%%
# BC case with helios TP
data_bc, bc_flux = read_in('BC', nsim, start_str = 'helios_tp_')

hcn_rain = [] # storing the rainout rates
for d in data_bc:
    hcn_rain.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))

rain = [] # storing the rainout rates
for d in data_bc:
    rain.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

#plot_rain(rain, bomb_rate, 'BC', figname = 'helios_tp_H2O_rainout_meteor.pdf', bc_flux_list = bc_flux, mol = 'Water', plot_Pearce = False, yscale = 'linear')
#plot_rain(hcn_rain, bomb_rate, 'BC', figname = 'helios_tp_HCN_rainout_meteor.pdf', bc_flux_list = bc_flux, plot_Pearce = False, yscale = 'linear')
#plot_vertical_n(data_bc, 'HCN', bomb_rate, 'BC', figname = 'helios_tp_HCN_air_meteor.pdf')
#plot_end_time(data_bc, figname = 'helios_tp_end_time_meteor.pdf')
#plot_evo_layer(data_bc, 'HCN', figname = 'helios_tp_hcn_evo_meteor.pdf')
#plot_convergence(data_bc, figname = 'helios_tp_convergence_meteor.pdf')
plot_rain_converged(data_bc, hcn_rain, bomb_rate, 'BC', figname = 'conv_BC_hcn_rain.pdf', rain_spec = 'HCN_rain', extra_list = bc_flux)
plot_rain_converged(data_bc, rain, bomb_rate, 'BC', figname = 'conv_BC_rain.pdf', rain_spec = 'H2O_rain', extra_list = bc_flux)
# %%
# C/O case

data_CtoO, C_to_O = read_in('CtoO', nsim)

hcn_rain_CtoO = []
for d in data_CtoO:
    hcn_rain_CtoO.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    
rain_CtoO = [] # storing the rainout rates
for d in data_CtoO:
    rain_CtoO.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
plot_rain(hcn_rain_CtoO, C_to_O, 'C_to_O', figname = 'HCN_rainout_C_to_O.pdf')
plot_rain(rain_CtoO, C_to_O, 'C_to_O', figname = 'H2O_rainout_C_to_O.pdf', mol = 'Water')
plot_vertical_n(data_CtoO, 'HCN', C_to_O, 'C_to_O', figname = 'HCN_air_C_to_O.pdf', f_size = 19, l_size = 18)
plot_end_time(data_CtoO, figname = 'end_time_C_to_O.pdf')
plot_evo_layer(data_CtoO, 'HCN', figname = 'hcn_evo_C_to_O.pdf')
plot_convergence(data_CtoO, figname = 'convergence_C_to_O.pdf')

#%%
# C/O case with HELIOS tP

data_CtoO, C_to_O = read_in('CtoO', nsim, start_str = 'helios_tp_')

hcn_rain_CtoO = []
for d in data_CtoO:
    hcn_rain_CtoO.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    
rain_CtoO = [] # storing the rainout rates
for d in data_CtoO:
    rain_CtoO.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
#plot_rain(hcn_rain_CtoO, C_to_O, 'C_to_O', figname = 'helios_tp_HCN_rainout_C_to_O.pdf', plot_Pearce = False)
#plot_rain(rain_CtoO, C_to_O, 'C_to_O', figname = 'helios_tp_H2O_rainout_C_to_O.pdf', mol = 'Water', plot_Pearce = False)
#plot_vertical_n(data_CtoO, 'HCN', C_to_O, 'C_to_O', figname = 'helios_tp_HCN_air_C_to_O.pdf', f_size = 19, l_size = 18)
#plot_end_time(data_CtoO, figname = 'helios_tp_end_time_C_to_O.pdf')
#plot_evo_layer(data_CtoO, 'HCN', figname = 'helios_tp_hcn_evo_C_to_O.pdf')
#plot_convergence(data_CtoO, figname = 'helios_tp_convergence_C_to_O.pdf')
plot_rain_converged(data_CtoO, hcn_rain_CtoO, C_to_O, 'CtoO', figname = 'conv_CtoO_hcn_rain.pdf', rain_spec = 'HCN_rain')
plot_rain_converged(data_CtoO, rain_CtoO, C_to_O, 'CtoO', figname = 'conv_CtoO_rain.pdf', rain_spec = 'H2O_rain')
plot_rain_converged(data_CtoO, hcn_rain_CtoO, C_to_O, 'CtoO', figname = 'conv_nonconv_CtoO_hcn_rain.pdf', rain_spec = 'HCN_rain', plot_non_conv = True)
plot_rain_converged(data_CtoO, rain_CtoO, C_to_O, 'CtoO', figname = 'conv_nonconv_CtoO_rain.pdf', rain_spec = 'H2O_rain', plot_non_conv = True)
# %%
for d in data_bc:
    print_max_re(d, 'H2')
# %%
# star case
import pandas as pd
data_star = read_in('star', number_of_sim = 13)
hcn_rain_star, rain_star = [], []
T_eff = pd.read_csv('/scratch/s2555875/stellar_flux/stellar_params.csv').T_eff

for d in data_star:
    hcn_rain_star.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    rain_star.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
#plot_rain(hcn_rain_star, T_eff, 'star', figname = 'HCN_rainout_star.pdf', plot_Pearce = False)
#plot_rain(rain_star, T_eff, 'star', figname = 'H2O_rainout_star.pdf', mol = 'Water', plot_Pearce = False)
#plot_vertical_n(data_star, 'HCN', T_eff, 'star', figname = 'HCN_air_star.pdf', f_size = 19, l_size = 18)
#plot_end_time(data_star, figname = 'end_time_star.pdf')
#plot_evo_layer(data_star, 'HCN', figname = 'hcn_evo_star.pdf')
#plot_convergence(data_star, figname = 'convergence_star.pdf')
plot_rain_converged(data_star, hcn_rain_star, T_eff, 'star', figname = 'conv_star_hcn_rain.pdf', rain_spec = 'HCN_rain')
plot_rain_converged(data_star, rain_star, T_eff, 'star', figname = 'conv_star_rain.pdf', rain_spec = 'H2O_rain')
plot_rain_converged(data_star, hcn_rain_star, T_eff, 'star', figname = 'conv_nonconv_star_hcn_rain.pdf', rain_spec = 'HCN_rain', plot_non_conv = True)
plot_rain_converged(data_star, rain_star, T_eff, 'star', figname = 'conv_nonconv_star_rain.pdf', rain_spec = 'H2O_rain', plot_non_conv = True)
# %%
# distance case
a_list = np.linspace(0.82, 1.4, nsim, endpoint = True) #HZ limits from Kopprapau et al. (2013) are 0.99 and 1.7, let's explore a bit more, from Venus to 2 au

data_dist, T_surf = read_in('dist', nsim)
hcn_rain_dist, rain_dist = [], []

for d in data_dist:
    hcn_rain_dist.append(rainout(d, rain_spec = 'HCN_rain', g_per_mol = 27))
    rain_dist.append(rainout(d, rain_spec = 'H2O_rain', g_per_mol = 18))

# do all the ploting
#plot_rain(hcn_rain_dist, a_list, 'dist', figname = 'HCN_rainout_dist.pdf', surf_temp = T_surf, plot_Pearce = False)
#plot_rain(rain_dist, a_list, 'dist', figname = 'H2O_rainout_dist.pdf', mol = 'Water', surf_temp = T_surf, plot_Pearce = False)
#plot_vertical_n(data_dist, 'HCN', a_list, 'dist', figname = 'HCN_air_dist.pdf', f_size = 19, l_size = 18)
#plot_end_time(data_dist, figname = 'end_time_dist.pdf')
#plot_evo_layer(data_dist, 'HCN', figname = 'hcn_evo_dist.pdf')
#plot_convergence(data_dist, figname = 'convergence_dist.pdf')
plot_rain_converged(data_dist, hcn_rain_dist, a_list, 'dist', figname = 'conv_dist_hcn_rain.pdf', rain_spec = 'HCN_rain', extra_list = T_surf)
plot_rain_converged(data_dist, rain_dist, a_list, 'dist', figname = 'conv_dist_rain.pdf', rain_spec = 'H2O_rain', extra_list = T_surf)
plot_rain_converged(data_dist, hcn_rain_dist, a_list, 'dist', figname = 'conv_nonconv_dist_hcn_rain.pdf', rain_spec = 'HCN_rain', extra_list = T_surf, plot_non_conv = True)
plot_rain_converged(data_dist, rain_dist, a_list, 'dist', figname = 'conv_nonconv_dist_rain.pdf', rain_spec = 'H2O_rain', extra_list = T_surf, plot_non_conv = True)
# %%
