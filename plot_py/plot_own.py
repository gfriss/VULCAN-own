#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import plot_reset as pr
import os
wd = os.getcwd()
os.chdir('../')
import parallel_functions as pf
os.chdir(wd)

pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6)

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
  yt = dat['variable']['y_time'][-800:,:,:]
  if years:
    tt = dat['variable']['t_time'][-800:] / (24*365.24*3600)
    ax[-1].set_xlabel('Time [yr]')
  else:
    tt = dat['variable']['t_time'][-800:]
    ax[-1].set_xlabel('Time [s]')
  
  dat_spec = dat['variable']['species']
  colour = ['r', 'orange', 'b', 'peru', 'c', 'k', 'magenta', 'purple', 'y', 'g', 'gray', 'yellow']
  lab = []
  if mol == None:
    mol = ['H2', 'N2', 'H2O', 'CO2', 'O2', 'CO', 'CH4', 'HCN', 'H2CO', 'OH', 'H2O_l_s', 'CH3']
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
vul_data_crahcno = scratch + 'output/archean.vul'
with open(vul_data_crahcno, 'rb') as handle:
  data_crahcno = pickle.load(handle)
  
vul_data_crahcno_ignore_5 = scratch + 'output/archean_ignore_5.vul'
with open(vul_data_crahcno_ignore_5, 'rb') as handle:
  data_crahcno_ignore_5 = pickle.load(handle)
  
vul_data_ncho = scratch + 'output/archean_ncho.vul'
with open(vul_data_ncho, 'rb') as handle:
  data_ncho = pickle.load(handle)
  
vul_data_ncho_ignore_5 = scratch + 'output/archean_ncho_ignore_5.vul'
with open(vul_data_ncho_ignore_5, 'rb') as handle:
  data_ncho_ignore_5 = pickle.load(handle)
  
vul_data_ncho_fix_bt_sp = scratch + 'output/archean_ncho_fix_bt_sp.vul'
with open(vul_data_ncho_fix_bt_sp, 'rb') as handle:
  data_ncho_fix_bt_sp = pickle.load(handle)
# %%
import plot_reset as pr
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6)
#%%
sp_to_plot = 'H2O_l_s'
fig, ax = plt.subplots(tight_layout = True)
#ax.plot(data_crahcno['variable']['ymix'][:,data_crahcno['variable']['species'].index(sp_to_plot)], data_crahcno['atm']['pco']/1.e6, label = 'CRAHCNO', color = 'blue', linestyle = '-')
#ax.plot(data_crahcno_ignore_5['variable']['ymix'][:,data_crahcno_ignore_5['variable']['species'].index(sp_to_plot)], data_crahcno_ignore_5['atm']['pco']/1.e6, label = 'CRAHCNO - ignore 5', linestyle = '--', color = 'orange')
ax.plot(data_ncho['variable']['ymix'][:,data_ncho['variable']['species'].index(sp_to_plot)], data_ncho['atm']['pco']/1.e6, label = 'NCHO', color = 'green', linestyle = '-')
#ax.plot(data_ncho_ignore_5['variable']['ymix'][:,data_ncho_ignore_5['variable']['species'].index(sp_to_plot)], data_ncho_ignore_5['atm']['pco']/1.e6, label = 'NCHO - ignore 5', linestyle = '--', color = 'red')
ax.plot(data_ncho_fix_bt_sp['variable']['ymix'][:,data_ncho_fix_bt_sp['variable']['species'].index(sp_to_plot)], data_ncho_fix_bt_sp['atm']['pco']/1.e6, label = 'fix bt sp', color = 'red', linestyle = '-')
ax.set_xlabel('X_{}'.format(sp_to_plot))
ax.set_ylabel('Pressure [bar]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((1e-12, None))
ax.invert_yaxis()
fig.legend(loc = (0.65,0.75))
#fig.savefig(scratch + 'plot/networks_{}.pdf'.format(sp_to_plot))
# %%
rain_spec = 'HCN_rain'
#dat_list = [data_crahcno, data_crahcno_ignore_5, data_ncho, data_ncho_ignore_5]
#name_list = ['CRAHCNO', 'CRAHCNO_ignore_5', 'NCHO', 'NCHO_ignore_5']
dat_list = [data_ncho, data_ncho_fix_bt_sp]
name_list = ['NCHO', 'NCHO_fix_bt_sp']
rain_rate = []
end_time = []
for dat in dat_list:
  rain_rate.append((np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'])*5.237e-13*(27/1000.)) # kg/m2yr
  end_time.append(dat['variable']['t'])
fig, ax = plt.subplots(tight_layout = True)
ax.scatter(name_list, rain_rate, c = 'red')
ax.set_ylabel('{} rate [kg/m2yr]'.format(rain_spec))
ax.set_yscale('log')
ax.set_ylim((min(rain_rate)/2, None))
#fig.savefig(scratch + 'plot/network_{}.pdf'.format(rain_spec))
# %%
fig, ax = plt.subplots(tight_layout = True)
ax.scatter(name_list, end_time, c = 'red')
ax.set_ylabel('End of sim time [s]')
ax.set_yscale('log')
ax.set_ylim((min(end_time)/2, None))
fig.savefig(scratch + 'plot/network_end_time.pdf')
# %%
import matplotlib.colors as mcolors
import matplotlib.cm as cm
pr.reset_plt(ticksize = 15, fontsize = 18, fxsize = 10, fysize = 8)
#%%
sp_to_plot = 'HCN'
fig, ax = plt.subplots(tight_layout = True)
markers = ['o', 's', 'v', '^']
norm = mcolors.LogNorm(vmin=min(data_crahcno['variable']['t_time']), vmax=max(data_crahcno['variable']['t_time']))
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)
for d,m,name in zip(dat_list,markers,name_list):
  t_time = d['variable']['t_time']
  for i in range(len(t_time)):
    if i == 0:
      ax.scatter(d['variable']['y_time'][i, ::5, d['variable']['species'].index(sp_to_plot)], d['atm']['pco'][::5]/1.e6, marker = m, label = name, c=mapper.to_rgba(t_time[i]*np.ones_like( d['atm']['pco'][::5]/1.e6)))
    elif i > 0 and i%20 == 0:
      ax.scatter(d['variable']['y_time'][i, ::5, d['variable']['species'].index(sp_to_plot)], d['atm']['pco'][::5]/1.e6, marker = m, c=mapper.to_rgba(t_time[i]*np.ones_like( d['atm']['pco'][::5]/1.e6)))
ax.set_xlabel('n_{} [cm^-3]'.format(sp_to_plot))
ax.set_ylabel('Pressure [bar]')
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlim((1e-2, 1e12))
#ax.set_xlim((None, 1e-21))
ax.invert_yaxis()
fig.legend(loc = (0.12,0.45))
fig.colorbar(mapper, label = 'Time [s]')
fig.savefig(scratch + 'plot/networks_evo_{}.pdf'.format(sp_to_plot))
# %%
# comparing synthetic and MUSCLES spectra
import numpy as np
import matplotlib.pyplot as plt

muscles = np.genfromtxt('/scratch/s2555875/stellar_flux/hd40307.txt', comments = '#', names = ['lambda', 'flux'])
#synthetic = np.genfromtxt('/scratch/s2555875/stellar_flux/sim_06_star.txt', comments = '#', names = ['lambda', 'flux'])
bosz_wl = np.genfromtxt('/scratch/s2555875/BOSZ_spectra/bosz2024_wave_r500.txt', comments = '#')
bosz_wl *= 1e-1 # nm
bosz_F = np.genfromtxt('/scratch/s2555875/BOSZ_spectra/bosz2024_mp_t5000_g+4.5_m+0.00_a+0.00_c+0.00_v0_r500_resam.txt', comments = '#', usecols = 0)
bosz_F *= 4*np.pi * 10 # erg/cm2/s/nm
nccs = np.genfromtxt('/scratch/s2555875/stellar_flux/hd22049.txt', comments = '#', names = ['lambda', 'flux'], skip_header = 16, skip_footer = 1)
plt.plot(muscles['lambda'], muscles['flux'], label = 'MUSCLES - HD40307')
plt.plot(nccs['lambda']/10, nccs['flux']*10*(0.53*1.5e11 / 5.1156e8)**2, label = 'NCCS - HD22049') # planet at 0.53 AU from star according to Seguro et al 2003
plt.plot(bosz_wl, bosz_F, label = 'BOSZ - 5000K')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Flux [ergs/cm^2/s/nm]')
plt.yscale('log')
plt.xscale('log')
plt.ylim((9e-2, 5e8))
plt.legend(loc = 'lower left')
plt.savefig('/scratch/s2555875/plot/syntheticnccs_muscles_T5000.pdf')
# %%
vul_data_ncho_1 = scratch + 'output/archean_ncho_1.vul'
with open(vul_data_ncho_1, 'rb') as handle:
  data_ncho_1 = pickle.load(handle)
  
vul_data_ncho_2 = scratch + 'output/archean_ncho_2.vul'
with open(vul_data_ncho_2, 'rb') as handle:
  data_ncho_2 = pickle.load(handle)
  
vul_data_ncho_3 = scratch + 'output/archean_ncho_3.vul'
with open(vul_data_ncho_3, 'rb') as handle:
  data_ncho_3 = pickle.load(handle)
  
vul_data_ncho_4 = scratch + 'output/archean_ncho_4.vul'
with open(vul_data_ncho_4, 'rb') as handle:
  data_ncho_4 = pickle.load(handle)
  
vul_data_ncho_5 = scratch + 'output/archean_ncho_5.vul'
with open(vul_data_ncho_5, 'rb') as handle:
  data_ncho_5 = pickle.load(handle)
  
vul_data_ncho_6 = scratch + 'output/archean_ncho_6.vul'
with open(vul_data_ncho_6, 'rb') as handle:
  data_ncho_6 = pickle.load(handle)
#%%
sp_to_plot = ['H2CO', 'CH3', 'CH4', 'H2O_l_s', 'HCN', 'H2', 'H', 'H2O', 'CO2', 'CO', 'OH']
for sp in sp_to_plot:
  fig, ax = plt.subplots(tight_layout = True)
  ax.plot(data_ncho_1['variable']['ymix'][:,data_ncho_1['variable']['species'].index(sp)], data_ncho_1['atm']['pco']/1.e6, label = '1', color = 'blue', linestyle = '-')
  ax.plot(data_ncho_2['variable']['ymix'][:,data_ncho_2['variable']['species'].index(sp)], data_ncho_2['atm']['pco']/1.e6, label = '2', linestyle = '--', color = 'orange')
  ax.plot(data_ncho_3['variable']['ymix'][:,data_ncho_3['variable']['species'].index(sp)], data_ncho_3['atm']['pco']/1.e6, label = '3', color = 'green', linestyle = '-')
  ax.plot(data_ncho_4['variable']['ymix'][:,data_ncho_4['variable']['species'].index(sp)], data_ncho_4['atm']['pco']/1.e6, label = '4', linestyle = '--', color = 'red')
  ax.plot(data_ncho_5['variable']['ymix'][:,data_ncho_5['variable']['species'].index(sp)], data_ncho_5['atm']['pco']/1.e6, label = '5', linestyle = ':', color = 'black')
  ax.plot(data_ncho_6['variable']['ymix'][:,data_ncho_6['variable']['species'].index(sp)], data_ncho_6['atm']['pco']/1.e6, label = '6', linestyle = '-.', color = 'magenta')
  ax.set_xlabel('X_{}'.format(sp))
  ax.set_ylabel('Pressure [bar]')
  ax.set_xscale('log')
  ax.set_yscale('log')
  #ax.set_xlim((1e-12, None))
  ax.invert_yaxis()
  fig.legend()
  fig.savefig(scratch + 'plot/ncho_rtoltest_{}.pdf'.format(sp))
#%%
rain_spec = ['H2O_rain', 'HCN_rain']
dat_list = [data_ncho_1, data_ncho_2, data_ncho_3, data_ncho_4, data_ncho_5, data_ncho_6]
name_list = ['1', '2', '3', '4', '5', '6']
for rain_sp in rain_spec:
  rain_rate = []
  end_time = []
  for dat in dat_list:
    rain_rate.append((np.sum(dat['variable']['y_rain'][rain_sp][:-1] * dat['atm']['dzi']) / dat['variable']['dt'])*5.237e-13*(27/1000.)) # kg/m2yr
    end_time.append(dat['variable']['t'])
  fig, ax = plt.subplots(tight_layout = True)
  ax.scatter(name_list, rain_rate, c = 'red')
  ax.set_ylabel('{} rate [kg/m2yr]'.format(rain_sp))
  ax.set_yscale('log')
  ax.set_ylim((min(rain_rate)/2, None))
  fig.savefig(scratch + 'plot/ncho_rtol_test_{}.pdf'.format(rain_sp))
#%%
# plot of TP and Eddy=diffusion profile for Archean
with open(scratch + 'output/archean_ncho.vul', 'rb') as handle:
  data_archean = pickle.load(handle)
colour = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig, ax = plt.subplots(tight_layout = True)
ax.plot(data_archean['atm']['Tco'], data_archean['atm']['pco']/1e6, label = 'T profile', color = colour[0], linestyle = '-')
ax.set_xlabel('Temperature [K]', color = colour[0])
ax.tick_params(axis = 'x', labelcolor = colour[0])
ax.set_ylabel('Pressure [bar]')
ax.set_yscale('log')
ax.invert_yaxis()
ax1 = ax.twiny()
ax1.plot(data_archean['atm']['Kzz'], data_archean['atm']['pco'][:-1]/1e6, label = 'Eddy diffusion', color = colour[6], linestyle = '--')
ax1.set_xlabel('Eddy diffusion [cm2/s]', color = colour[6])
ax1.set_xscale('log')
ax1.tick_params(axis = 'x', labelcolor = colour[6])
fig.savefig(scratch + 'plot/TPs/archean_T_Kzz.pdf')
#%%
original_archean = scratch + 'output/archean_ncho_nowash.vul'
with open(original_archean, 'rb') as handle:
  data_original_archean = pickle.load(handle)
  
original_archean_updated = scratch + 'output/archean_updated_ncho_nowash.vul'
with open(original_archean_updated, 'rb') as handle:
  data_original_archean_updated = pickle.load(handle)  
  
stellarcut_archean = scratch + 'output/archean_ncho_stellarcut.vul'
with open(stellarcut_archean, 'rb') as handle:
  data_stellarcut_archean = pickle.load(handle)

stellarcut_archean_updated = scratch + 'output/archean_updated_ncho_stellarcut.vul'
with open(stellarcut_archean_updated, 'rb') as handle:
  data_stellarcut_archean_updated = pickle.load(handle)
  
sims_to_compare = ['Archean', 'Archean updated', 'Archean stellar cut', 'Archean updated stellar cut']
sim_data = [data_original_archean, data_original_archean_updated, data_stellarcut_archean, data_stellarcut_archean_updated]
  
pr.reset_plt(ticksize = 20, fontsize = 24, fxsize = 16, fysize = 8)
  
fig, ax = plt.subplots(tight_layout = True, ncols = 2, sharey = True, sharex = True)
ax = ax.flatten()
for i, sim in enumerate(sim_data):
  ax[0].plot(sim['variable']['ymix'][:,sim['variable']['species'].index('HCN')], sim['atm']['pco']/1.e6, label = sims_to_compare[i])
  ax[1].plot(sim['variable']['ymix'][:,sim['variable']['species'].index('H2O_l_s')], sim['atm']['pco']/1.e6, label = sims_to_compare[i])
ax[0].set_xlabel(r'$X_{HCN}$')
ax[0].set_ylabel('Pressure [bar]')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlim((1e-21, None))
ax[0].invert_yaxis()
ax[1].set_xlabel(r'$X_{H_2O_{l/s}}$')
ax[1].legend()
fig.savefig(scratch + 'plot/vertical_profiles/archean_stellarcut_compare.pdf')
#%%
pr.reset_plt(ticksize = 13, fontsize = 15, fxsize = 8, fysize = 6)

hcn_rain, h2o_rain = [], []
for d in [data_original_archean, data_original_archean_updated, data_stellarcut_archean, data_stellarcut_archean_updated]:
  hcn_rain.append(pf.rainout(dat = d, rain_spec = 'HCN_rain', g_per_mol = 27))
  h2o_rain.append(pf.rainout(dat = d, rain_spec = 'H2O_rain', g_per_mol = 18))
fig, ax = plt.subplots(tight_layout = True)
ax.plot(sims_to_compare, hcn_rain, linestyle = '', marker = 'o')
ax.set_yscale('log')
ax.tick_params(axis='x', rotation=15)
ax.set_ylabel('HCN rainout [kg/m2/yr]')
ax2 = ax.twinx()
ax2.plot(sims_to_compare, h2o_rain, linestyle = '', marker = 's', color = 'orange', alpha = 0.7)
ax2.set_yscale('log')
ax2.set_ylabel('H2O rainout [kg/m2/yr]')
fig.savefig(scratch + 'plot/rainout_rates/archean_rainout_compare.pdf')
# %%