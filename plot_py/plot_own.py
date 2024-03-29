#%%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pickle

datastore = '/tmp/datastore/s2555875/'
vul_own = 'VULCAN-own/'
vul_mas = 'VULCAN-master/'

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

vul_crahcno_ox = datastore + vul_own + 'output/CRAHCNO_ox.vul'
with open(vul_crahcno_ox, 'rb') as handle:
  data_CRAHCNO_ox = pickle.load(handle)

vul_crahcno_red = datastore + vul_mas + 'output/CRAHCNO_red.vul'
with open(vul_crahcno_red, 'rb') as handle:
  data_CRAHCNO_red = pickle.load(handle)

vulcan_spec = data_CRAHCNO_ox['variable']['species']
spec = ['HCN', 'H2O_l_s']
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
def plot_evo_layer(yt, tt, n, mol, xscale = 'log', yscale = 'log', ylim = None, figname = None):
  fig, ax = plt.subplots()
  ytot = np.sum(yt[:, n, data_CRAHCNO_ox['atm']['gas_indx']], axis = 1)
  for molec in mol:
    ax.plot(tt/(24*365.24*3600), yt[:, n, vulcan_spec.index(molec)]/ytot, label = molec)
  ax.set_xlabel('Time [yr]')
  ax.set_ylabel('Mixing ratio')#(r'n [cm$^{-3}$]')
  ax.set_yscale(yscale)
  ax.set_xscale(xscale)
  ax.legend()
  if ylim != None:
    ax.set_ylim(ylim)
  if figname != None:
    fig.savefig(figname + '.png')
# %%
