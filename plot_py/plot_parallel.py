#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle

nsim = 15 # hardcoded for now, change later...
bomb_rate = np.linspace(3e23, 1e25, nsim) # values from Pearce et al. (2022) Fig A4 min and max
out_folder = '/scratch/s2555875/output/'
bc_folder= '/scratch/s2555875/BC_files/'
plot_folder = '/scratch/s2555875/plot/'
data = [] # list to store the results (dicts)
bc_spec = ['H2']#, 'CO2', 'CH4']
bc_flux = {}
for sp in bc_spec:
    bc_flux[sp] = []
bc_linestyle = ['-']#, '--', ':']

for i in range(nsim):
    sim = 'sim_' # setting up the names of files
    if i < 10:
        sim += '0' + str(i)
    else:
        sim += str(i)
    sim += '_onlyH2.vul'
    
    with open(out_folder+sim, 'rb') as handle: # reading in files
        data.append(pickle.load(handle))

    with open(bc_folder+'BC_bot_'+sim[:-11]+'.txt') as f: # vas sim[:-3]+'txt
        count = 0
        for line in f:
            lin = line.split()
            if line[0] != '#' and lin[0] in bc_spec:
                bc_flux[lin[0]].append(float(lin[1]))
                count +=1
            if count == len(bc_spec):
                break

# post simulations...
nsim_extra = 8
data_extra = []
bc_flux_extra = {}
for sp in bc_spec:
    bc_flux_extra[sp] = []
for i in range(nsim_extra):
    sim = 'sim_' # setting up the names of files
    if i < 10:
        sim += '00' + str(i)
    else:
        sim += str(i)
    sim += '_onlyH2.vul'
    
    with open(out_folder+sim, 'rb') as handle: # reading in files
        data_extra.append(pickle.load(handle))

    with open(bc_folder+'BC_bot_'+sim[:-11]+'.txt') as f: # vas sim[:-3]+'txt
        count = 0
        for line in f:
            lin = line.split()
            if line[0] != '#' and lin[0] in bc_spec:
                bc_flux_extra[lin[0]].append(float(lin[1]))
                count +=1
            if count == len(bc_spec):
                break

data = [data[0]] + data_extra + data[1:] # to have bombardment rates in order
bc_flux['H2'] = [bc_flux['H2'][0]] + bc_flux_extra['H2'] + bc_flux['H2'][1:]
bomb_rate_extra = np.linspace(3.5e23, 9e23, nsim_extra)
bomb_rate = np.concatenate((np.concatenate((np.array([bomb_rate[0]]), bomb_rate_extra)), bomb_rate[1:]))
vulcan_spec = data[0]['variable']['species'] # globally st up the species (same network so can do)

def rainout(dat, rain_spec = 'HCN_rain', g_per_mol = 27):
    rain_rate = np.sum(dat['variable']['y'][:-1,vulcan_spec.index(rain_spec)] * dat['atm']['dzi']) # 1/cm2s
    rain_rate = rain_rate * 2.259e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

hcn_rain = [] # storing the rainout rates
for d in data:
    hcn_rain.append(rainout(d, rain_spec = 'H2O_l_s', g_per_mol = 18))

# plot rainout rates and BC conditions as a function of meteoritic mass delivery rate
# plus include results from Pearce at al. (2022)

from matplotlib.lines import Line2D

def create_dummy_line(**kwds):
    return Line2D([], [], **kwds)

with open(out_folder+'B_nofix.vul', 'rb') as handle:
    data_B = pickle.load(handle)

hcn_rain_B = 3.4e-12 #rainout(data_B)
bc_B = [2.3e+10, 3.0e11, 6.8e8]
bomb_B = 1.2e24 * 2.3/3.42 # maybe it is not scaled or my calc is bad
mark = ['*', 's', 'p']
#%%
fig, ax = plt.subplots(tight_layout = True)
colour = 'tab:blue'
ax.plot(bomb_rate, hcn_rain, color = colour, label = 'Own')
#ax.plot(bomb_B, hcn_rain_B, marker = '*', linestyle = '', label = 'Pearce et al. (2022)')
ax.set_xlabel(r'Mass delivery rate [g Gyr$^{-1}$]', fontsize=13)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'HCN rain-out rate [kg m$^{-2}$ yr$^{-1}$]', color = colour, fontsize=13)
ax.tick_params(axis = 'y', labelcolor = colour, labelsize = 12)
ax.legend(loc = 'center left')
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(12) 
ax1 = ax.twinx()
colour = 'tab:orange'
ax1.tick_params(axis = 'y', labelcolor = colour)
for i in range(len(bc_spec)):
    ax1.plot(bomb_rate, bc_flux[bc_spec[i]], color = colour, linestyle = bc_linestyle[i])
    #ax1.plot(bomb_B, bc_B[i], marker = mark[i], color = colour, linestyle = bc_linestyle[i], label = bc_spec[i])
ax1.set_ylabel(r'H2 flux from surface [cm$^{-2}$ s$^{-1}$]', color = colour, fontsize=13)
ax1.set_yscale('log')
#ax1.legend(loc = 'center right')
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(12) 

ax.tick_params(which='both', direction='out', width=1, length = 4)
ax1.tick_params(which='both', direction='out', width=1, length = 4)
ax1.tick_params(axis='both',labelsize=12)

#fig.savefig(plot_folder + 'HCN_rainout_meteor_onlyH2.pdf')
fig.savefig(plot_folder + 'H2O_rainout_meteor_onlyH2.png')

fig, ax = plt.subplots()
for i in range(len(data)):
    ax.plot(data[i]['variable']['y'][:, data[i]['variable']['species'].index('HCN')], data[i]['atm']['zco'][1:]/1e5, label = r'M${_del}% = ' + '{:.2e} g/Gyr'.format(bomb_rate[i]))
ax.set_xscale('log')
ax.set_xlabel(r'n [cm$^{-3}$]')
ax.set_ylabel('Height [km]')
ax.legend()

fig.savefig(plot_folder + 'HCN_air_meteoritic.png')
#%%
for i in range(len(data)):
    plt.plot(i, data[i]['variable']['t'], linestyle = '', marker = 'o', color = 'red')
plt.yscale('log')
plt.yscale('log')
plt.xlabel('Simulation number')
plt.ylabel('End of simulation time [s]')
plt.savefig(plot_folder + 'end_time_meteor.pdf')
#%%
for i in range(len(data)):
    plt.plot(i, data[i]['variable']['y'][0, data[i]['variable']['species'].index('HCN')], linestyle = '', marker = 'o')
plt.yscale('log')
plt.xlabel('Simulation number')
plt.ylabel('HCN number density [cm-3]')
#%%
plt.figure(figsize = (12, 8))
for i in range(len(data)):
    plt.plot(data[i]['variable']['t_time'], data[i]['variable']['y_time'][:, 0, data[i]['variable']['species'].index('HCN')], label = i)
plt.xlabel('Time [s]')
plt.ylabel('n [cm-3]')
plt.yscale('log')
plt.legend()
plt.ylim(1e-2,None)
plt.savefig(plot_folder + 'hcn_evo_meteor.pdf')
#%%
yconv_cri = 0.01
yconv_min = 0.1
slope_cri = 1.e-4
#plt.figure(figsize = (12,8))
for i in range(len(data)):
    longdy = data[i]['variable']['longdy']
    longdydt = data[i]['variable']['longdydt']
    slope_min = min( np.amin(data[i]['atm']['Kzz']/(0.1*data[i]['atm']['Hp'][:-1])**2) , 1.e-8)
    slope_min = max(slope_min, 1.e-10)
    if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min):
        plt.plot(i, 'Yes', 'ro')
    else:
        plt.plot(i, 'No', 'ro')
plt.xlabel('Simulation number')
plt.savefig(plot_folder + 'convergence_meteor.pdf')
#%%
def calc_rate(dat, spec_list, re_id, n):
    rate = dat['variable']['k'][re_id][n]
    for sp in spec_list:
        if sp != 'M':
            rate *= dat['variable']['y'][n, dat['variable']['species'].index(sp)]
    return rate

def plot_max_re(dat, mol, n = 0, prod = False):
    reaction = ''
    reagents_spec, products_spec = [], []
    k_rea = 0
    re_rate = 0
    for k,v in dat['variable']['Rf'].items():
        reagents_products = v.split('->') # separating production and destruction
        reagents = reagents_products[0]
        reagents_spec = reagents.split('+')
        reagents_spec = [r.strip() for r in reagents_spec] # stripping them from white spaces
        products = reagents_products[1].strip()
        products_spec = products.split('+')
        products_spec = [p.strip() for p in products_spec] # stripping them from white spaces
        if prod == False and mol in reagents_spec: # destruction so it is on the reagents side
            re_rate = calc_rate(dat, reagents_spec, k, n)
        elif prod == True and mol in products_spec: # production side
            re_rate = calc_rate(dat, products_spec, k, n)
        if re_rate > k_rea:
            k_rea = re_rate
            reaction = v
    print('The highest reaction rate in layer {} is k = {:.3e} cm-3s-1 for the reaction {}.'.format(n, k_rea, reaction))
# %%
