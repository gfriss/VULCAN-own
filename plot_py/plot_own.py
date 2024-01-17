#%%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pickle

vul_data = '../output/Earth.vul'

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

vulcan_spec = data['variable']['species']
spec = ['HCN', 'H2O_l_s']

vul_data_rainout = '../output/separate.vul'

with open(vul_data_rainout, 'rb') as handle:
  data_rainout = pickle.load(handle)

vul_data_h2o_bot = '../output/Earth_h2o_bot.vul'
with open(vul_data_h2o_bot, 'rb') as handle:
  data_h2o_bot = pickle.load(handle)

vul_data_frac_to_rate = '../output/frac_to_rate.vul'
with open(vul_data_frac_to_rate, 'rb') as handle:
  data_frac_to_rate = pickle.load(handle)

vul_data_low_henrys = '../output/low_henrys.vul'
with open(vul_data_low_henrys, 'rb') as handle:
  data_low_henrys = pickle.load(handle)

vul_data_rain_from_l_s = '../output/rain_from_l_s.vul'
with open(vul_data_rain_from_l_s, 'rb') as handle:
  data_rain_from_l_s = pickle.load(handle)

vul_data_corrected_list_use = '../output/corrected_list_use.vul'
with open(vul_data_corrected_list_use, 'rb') as handle:
  data_corrected_list_use = pickle.load(handle)

vul_data_corrected_high_henrys = '../output/corrected_high_henrys.vul'
with open(vul_data_corrected_high_henrys, 'rb') as handle:
  data_corrected_high_henrys = pickle.load(handle)

vul_data_corrected_mod = '/home/s2555875/VULCAN-master/output/mod_rate_calc.vul'
with open(vul_data_corrected_mod, 'rb') as handle:
  data_corrected_mod = pickle.load(handle)
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
plt.savefig('../plot/corrected_h2o_hcn_mixing.png')
#%%
plt.plot(data_corrected_list_use['variable']['y'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['y'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 10$')
plt.plot(data_corrected_high_henrys['variable']['y'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['y'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 100000$')

plt.xlabel(r'$n_{rainout}-n$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/compare_henrys_const_number_dens.png')
#%%
plt.plot(data_corrected_list_use['variable']['ymix'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 10$')
plt.plot(data_corrected_high_henrys['variable']['ymix'][:,vulcan_spec.index('HCN')]-data_h2o_bot['variable']['ymix'][:,vulcan_spec.index('HCN')], data_h2o_bot['atm']['zco'][1:]/1.e5, label = r'$KH = 100000$')

plt.xlabel(r'$X_{rainout}-X$')
plt.ylabel('Height [km]')
plt.legend()
plt.savefig('../plot/compare_henrys_const_mixing.png')
# %%
