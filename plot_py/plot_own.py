#%%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pickle

datastore = '/tmp/datastore/s2555875/'
vul_own = 'VULCAN-own/'
vul_mas = 'VULCAN-master/'
scratch = '/scratch/s2555875/'
#%%

vul_data = datastore + vul_own + 'output/Earth.vul'

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

vul_data_rainout = datastore + vul_own + 'output/separate.vul'

with open(vul_data_rainout, 'rb') as handle:
  data_rainout = pickle.load(handle)

vul_data_h2o_bot = datastore + vul_own + 'output/Earth_h2o_bot.vul'
with open(vul_data_h2o_bot, 'rb') as handle:
  data_h2o_bot = pickle.load(handle)

vul_data_frac_to_rate = datastore + vul_own + 'output/frac_to_rate.vul'
with open(vul_data_frac_to_rate, 'rb') as handle:
  data_frac_to_rate = pickle.load(handle)

vul_data_low_henrys = datastore + vul_own + 'output/low_henrys.vul'
with open(vul_data_low_henrys, 'rb') as handle:
  data_low_henrys = pickle.load(handle)

vul_data_rain_from_l_s = datastore + vul_own + 'output/rain_from_l_s.vul'
with open(vul_data_rain_from_l_s, 'rb') as handle:
  data_rain_from_l_s = pickle.load(handle)

vul_data_corrected_list_use = datastore + vul_own + 'output/corrected_list_use.vul'
with open(vul_data_corrected_list_use, 'rb') as handle:
  data_corrected_list_use = pickle.load(handle)

vul_data_corrected_high_henrys = datastore + vul_own + 'output/corrected_high_henrys.vul'
with open(vul_data_corrected_high_henrys, 'rb') as handle:
  data_corrected_high_henrys = pickle.load(handle)

vul_data_corrected_mod = datastore + vul_mas + 'output/mod_rate_calc.vul'
with open(vul_data_corrected_mod, 'rb') as handle:
  data_corrected_mod = pickle.load(handle)

vul_data_cor_mod_top = datastore + vul_mas + 'output/mod_with_top.vul'
with open(vul_data_cor_mod_top, 'rb') as handle:
  data_cor_mod_top = pickle.load(handle)

vul_data_cor_mod_bot = datastore + vul_mas + 'output/mod_with_bot.vul'
with open(vul_data_cor_mod_bot, 'rb') as handle:
  data_cor_mod_bot = pickle.load(handle)

vul_data_washout = datastore + vul_own + 'output/washout.vul'
with open(vul_data_washout, 'rb') as handle:
  data_washout = pickle.load(handle)

vul_data_Pearce_A_like = datastore + vul_own + 'output/Pearce_A_like.vul'
with open(vul_data_Pearce_A_like, 'rb') as handle:
  data_Pearce_A_like = pickle.load(handle)

vul_data_refined_param = datastore + vul_own + 'output/refined_param.vul'
with open(vul_data_refined_param, 'rb') as handle:
  data_refined_param = pickle.load(handle)

vul_data_refined_washout_raintrack = datastore + vul_own + 'output/refined_washout_raintrack.vul'
with open(vul_data_refined_washout_raintrack, 'rb') as handle:
  data_refined_washout_raintrack = pickle.load(handle)

vul_data_Pearce_oxidising = datastore + vul_own + 'output/Pearce_oxidising.vul'
with open(vul_data_Pearce_oxidising, 'rb') as handle:
  data_Pearce_oxidising = pickle.load(handle)

vul_data_Pearce_oxidising_mixtable = datastore + vul_own + 'output/Pearce_oxidising_mixtable.vul'
with open(vul_data_Pearce_oxidising_mixtable, 'rb') as handle:
  data_Pearce_oxidising_mixtable = pickle.load(handle)

vul_data_Pearce_oxidising_mixtable_no_settling = datastore + vul_own + 'output/Pearce_oxidising_mixtable_no_settling.vul'
with open(vul_data_Pearce_oxidising_mixtable_no_settling, 'rb') as handle:
  data_Pearce_oxidising_mixtable_no_settling = pickle.load(handle)

vul_data_H2_escape_Pearce_A = datastore + vul_mas + 'output/H2_escape_Pearce_A.vul'
with open(vul_data_H2_escape_Pearce_A, 'rb') as handle:
  data_H2_escape_Pearce_A = pickle.load(handle)
#%%
vul_crahcno_ox = datastore + vul_own + 'output/CRAHCNO_ox.vul'
with open(vul_crahcno_ox, 'rb') as handle:
  data_CRAHCNO_ox = pickle.load(handle)

vul_crahcno_red = datastore + vul_mas + 'output/CRAHCNO_red.vul'
with open(vul_crahcno_red, 'rb') as handle:
  data_CRAHCNO_red = pickle.load(handle)

vul_longrun_B = datastore + vul_own + 'output/longrun_B.vul'
with open(vul_longrun_B, 'rb') as handle:
  data_longrun_B = pickle.load(handle)

vul_longrun_B_Earthlike_Kzz = datastore + vul_own + 'output/longrun_B_Earthlike_Kzz.vul'
with open(vul_longrun_B_Earthlike_Kzz, 'rb') as handle:
  data_longrun_B_Earthlike_Kzz = pickle.load(handle)

#%%
vul_B_Pearce = scratch + 'output/B_Pearce.vul'
with open(vul_B_Pearce, 'rb') as handle:
  data_B_Pearce = pickle.load(handle)

vul_B_Pearce_long = scratch + 'output/B_Pearce_long.vul'
with open(vul_B_Pearce_long, 'rb') as handle:
  data_B_Pearce_long = pickle.load(handle)

vul_B_nofix = scratch + 'output/B_nofix.vul'
with open(vul_B_nofix, 'rb') as handle:
  data_B_nofix = pickle.load(handle)

vul_B_nofix_long = scratch + 'output/B_nofix_long.vul'
with open(vul_B_nofix_long, 'rb') as handle:
  data_B_nofix_long = pickle.load(handle)

vul_B_nofix_norain = scratch + 'output/B_nofix_norain.vul'
with open(vul_B_nofix_norain, 'rb') as handle:
  data_B_nofix_norain = pickle.load(handle)

vulcan_spec = data_B_Pearce['variable']['species']
spec = ['HCN', 'H2O_l_s']
#%%
vul_Earth_nofix_norain = scratch + 'output/Earth_nofix_norain.vul'
with open(vul_Earth_nofix_norain, 'rb') as handle:
  data_Earth_nofix_norain = pickle.load(handle)

vul_Earth_nofix_with_rain = scratch + 'output/Earth_nofix_with_rain.vul'
with open(vul_Earth_nofix_with_rain, 'rb') as handle:
  data_Earth_nofix_with_rain = pickle.load(handle)

vul_Earth_nofix_high_henries = scratch + 'output/Earth_nofix_high_henries.vul'
with open(vul_Earth_nofix_high_henries, 'rb') as handle:
  data_Earth_nofix_high_henries = pickle.load(handle)

vul_Earth_nofix_medium_henries = scratch + 'output/Earth_nofix_medium_henries.vul'
with open(vul_Earth_nofix_medium_henries, 'rb') as handle:
  data_Earth_nofix_medium_henries = pickle.load(handle)

vulcan_spec = data_Earth_nofix_norain['variable']['species']
spec = ['HCN', 'H2O_l_s']
#%%
fig, ax = plt.subplots(tight_layout = True, figsize = (10,6))
ax.plot(data_Earth_nofix_with_rain['variable']['y'][:,vulcan_spec.index('HCN')]-data_Earth_nofix_norain['variable']['y'][:,vulcan_spec.index('HCN')], data_Earth_nofix_with_rain['atm']['zco'][1:]/1.e5, label = r'$H = 10$')
ax.plot(data_Earth_nofix_high_henries['variable']['y'][:,vulcan_spec.index('HCN')]-data_Earth_nofix_norain['variable']['y'][:,vulcan_spec.index('HCN')], data_Earth_nofix_high_henries['atm']['zco'][1:]/1.e5, label = r'$H = 100000$')
ax.set_xlabel(r'$n_{rainout}-n_0$', fontsize = 18)
#ax.set_xscale('symlog')
ax.set_ylabel('Height [km]', fontsize = 18)
ax.tick_params(axis = 'y', labelsize = 17)
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(17) 
ax.legend(fontsize = '16')
ax1 = ax.twinx()
ax1.set_ylabel('Pressure [bar]', fontsize = 18)
ax1.set_ylim(2, 1e-8)
ax1.set_yscale('log')
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(17) 

ax.tick_params(which='both', direction='out', width=1, length = 4)
ax1.tick_params(which='both', direction='out', width=1, length = 4)
ax1.tick_params(axis='y',labelsize=17)
ax.tick_params(axis='x',labelsize=17)
fig.savefig(scratch+'plot/hcn_rain_chem.pdf')
#%%
for sp in spec:
  plt.plot(data_rainout['variable']['ymix'][:,vulcan_spec.index(sp)]-data_rainout['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/separate.pdf')

#%%
for sp in spec:
  plt.plot(data_frac_to_rate['variable']['ymix'][:,vulcan_spec.index(sp)]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
#plt.savefig('../plot/separate.pdf')
#%%
for sp in spec:
  plt.plot(data_rain_from_l_s['variable']['ymix'][:,vulcan_spec.index(sp)]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
# %%
vul_data_no_botflux = '../output/no_h2o_bot.vul'

with open(vul_data_no_botflux, 'rb') as handle:
  data_no_botflux = pickle.load(handle)
#%%
for sp in spec:
  plt.plot(data_no_botflux['variable']['ymix'][:,vulcan_spec.index(sp)]-data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/separate_no_botflux.pdf')
# %%
vul_data_no_botflux_hh = '../output/no_h2o_bot_high_henrys.vul'

with open(vul_data_no_botflux_hh, 'rb') as handle:
  data_no_botflux_hh = pickle.load(handle)
#%% # did not help...
for sp in spec:
  plt.plot(data_no_botflux_hh['variable']['ymix'][:,vulcan_spec.index(sp)]-data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()

# %%
for sp in spec:
  plt.plot(data_corrected_list_use['variable']['ymix'][:,vulcan_spec.index(sp)]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/corrected_list_use.pdf')

#%%
for sp in spec:
  plt.plot(data_corrected_high_henrys['variable']['ymix'][:,vulcan_spec.index(sp)]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index(sp)], data_h2o_bot['atm']['zco'][1:]/1.e5, label = sp)

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/corrected_h2o_hcn_mixing.pdf')
#%%
plt.plot(data_corrected_list_use['variable']['y'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['y'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 10$')
plt.plot(data_corrected_high_henrys['variable']['y'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['y'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 100000$')

plt.xlabel(r'$n_{rainout}-n$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/compare_henrys_const_number_dens.pdf')
#%%
plt.plot(data_corrected_list_use['variable']['ymix'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 10$')
plt.plot(data_corrected_high_henrys['variable']['ymix'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 100000$')

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/compare_henrys_const_mixing.pdf')
# %%
fig, ax1 = plt.subplots()

ax1.plot(data_corrected_list_use['variable']['y'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = 'No washout')
ax1.set_xlabel(r'$n_{scavenging}$')
ax1.set_ylabel('Height [km]')

ax2 = ax1.twiny()
ax2.plot(data_corrected_list_use['variable']['y'][:,vulcan_spec.index('HCN')]-data_washout['variable']['y'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, color = 'orange', linestyle = '--', label = 'Difference')
ax2.set_xlabel(r'$n_{scavenging} - n_{washout}$')
fig.legend(loc = (0.71,0.76))
fig.savefig('../plot/hcn_washout_compare.pdf')
#%%
def plot_time_evo(yt, tt, n, mol, xscale = 'log', yscale = 'log', ylim = None):
  alpha = np.linspace(0.1, 1, n)

  fig, ax = plt.subplots()
  for layer in range(n):
    if layer in [0, n//4, n//2, 3*n//4, n-1]:
      ax.plot(tt, yt[:, layer, vulcan_spec.index(mol)], label = 'layer = {}'.format(layer), alpha = alpha[layer], color = 'red')
    else:
      ax.plot(tt, yt[:, layer, vulcan_spec.index(mol)], alpha = alpha[layer], color = 'red')

  ax.set_xlabel('Time [s]')
  ax.set_ylabel(r'n [cm$^{-3}$]')
  ax.set_yscale(yscale)
  ax.set_xscale(xscale)
  ax.legend()
  if ylim != None:
    ax.set_ylim(ylim)
# %%
def plot_evo_layer(dat, n, mol = None, xscale = 'log', yscale = 'log', ylim = None, figname = None, years = False):
  fig, ax = plt.subplots()
  yt = dat['variable']['y_time']
  if years:
    tt = dat['variable']['t_time'] / (24*365.24*3600)
    ax.set_xlabel('Time [yr]')
  else:
    tt = dat['variable']['t_time']
    ax.set_xlabel('Time [s]')
  
  dat_spec = dat['variable']['species']
  ytot = np.sum(yt[:, n, dat['atm']['gas_indx']], axis = 1)
  colour = ['r', 'orange', 'b', 'peru', 'c', 'k', 'magenta', 'purple', 'y', 'g']
  if mol == None:
    mol = ['H2', 'N2', 'H2O', 'CO2', 'O2', 'CO', 'CH4', 'HCN', 'H2CO', 'OH']

  for molec,c in zip(mol, colour):
    ax.plot(tt, yt[:, n, dat_spec.index(molec)]/ytot, label = molec, color = c)
  ax.set_ylabel('Mixing ratio')#(r'n [cm$^{-3}$]')
  ax.set_yscale(yscale)
  ax.set_xscale(xscale)
  ax.legend()
  if ylim != None:
    ax.set_ylim(ylim)
  if figname != None:
    fig.savefig(scratch + 'plot/' + figname + '.png')

def plot_many_layers(dat, n, mol = None, xscale = 'log', yscale = 'log', ylim = None, figname = None, years = False):
  fig, ax = plt.subplots(nrows = len(n), sharex = True, sharey = True, figsize = (8,4*len(n)), tight_layout = True)
  yt = dat['variable']['y_time']
  if years:
    tt = dat['variable']['t_time'] / (24*365.24*3600)
    ax[-1].set_xlabel('Time [yr]')
  else:
    tt = dat['variable']['t_time']
    ax[-1].set_xlabel('Time [s]')
  
  dat_spec = dat['variable']['species']
  colour = ['r', 'orange', 'b', 'peru', 'c', 'k', 'magenta', 'purple', 'y', 'g']
  lab = []
  if mol == None:
    mol = ['H2', 'N2', 'H2O', 'CO2', 'O2', 'CO', 'CH4', 'HCN', 'H2CO', 'OH']
  for i in range(len(ax)):
    ytot = np.sum(yt[:, n[-(i+1)], dat['atm']['gas_indx']], axis = 1)
    for molec,c in zip(mol, colour):
      if i == len(ax) - 1:
        ax[-1].plot(tt, yt[:, n[-(i+1)], dat_spec.index(molec)]/ytot, label = molec, color = c)
        lab.append(molec)
      else:
        ax[i].plot(tt, yt[:, n[-(i+1)], dat_spec.index(molec)]/ytot, color = c)
    ax[i].set_ylabel('Mixing ratio')#(r'n [cm$^{-3}$]')
    ax[i].set_yscale(yscale)
    ax[i].set_xscale(xscale)
    ax[i].set_title('Layer ' + str(n[-(i+1)]))
    #ax[i].set_xlim((1e12, None))
  if ylim != None:
      ax[0].set_ylim(ylim)
  fig.legend(labels = lab, bbox_to_anchor = (1.12, 0.65))
  if figname != None:
    fig.savefig(scratch + 'plot/' + figname + '.png', bbox_inches = 'tight')
# %%
# calc rain out rate # below gives 1.23e-8 which is much higher than for Ben...
hcn_rain_rate = np.sum(data_B_nofix['variable']['y_rain']['HCN_rain'][:-1] * data_B_nofix['atm']['dzi']) / data_B_nofix['variable']['dt'] # 1/cm2s
print(hcn_rain_rate * 1.24e-14) # kg/m2yr
#%%
# rain similar to Pearce et al (2022) without my rain, using the deposition velocity
hcn_rain_rate = np.sum(data_B_nofix_norain['variable']['y'][:,vulcan_spec.index('HCN')]) * data_B_nofix_norain['atm']['bot_vdep'][vulcan_spec.index('HCN')] # 1/cm2s
print(hcn_rain_rate * 1.24e-14) # kg/m2yr
# %%
