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

def rainout(dat):
    hcn_rain_rate = np.sum(dat['variable']['y'][:-1,vulcan_spec.index('HCN_rain')] * dat['atm']['dzi']) # 1/cm2s
    return hcn_rain_rate * 1.24e-14 # kg/m2yr

hcn_rain = [] # storing the rainout rates
for d in data:
    hcn_rain.append(rainout(d))

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

fig, ax = plt.subplots(tight_layout = True)
colour = 'tab:blue'
ax.plot(bomb_rate, hcn_rain, color = colour, label = 'Own')
ax.plot(bomb_B, hcn_rain_B, marker = '*', linestyle = '', label = 'Pearce et al. (2022)')
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
    ax1.plot(bomb_B, bc_B[i], marker = mark[i], color = colour, linestyle = bc_linestyle[i], label = bc_spec[i])
ax1.set_ylabel(r'H2 flux from surface [cm$^{-2}$ s$^{-1}$]', color = colour, fontsize=13)
ax1.set_yscale('log')
#ax1.legend(loc = 'center right')
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(12) 

ax.tick_params(which='both', direction='out', width=1, length = 4)
ax1.tick_params(which='both', direction='out', width=1, length = 4)
ax1.tick_params(axis='both',labelsize=12)

fig.savefig(plot_folder + 'HCN_rainout_meteor_onlyH2.pdf')

