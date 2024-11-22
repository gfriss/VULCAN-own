#%%
import numpy as np
import matplotlib.pyplot as plt
#import plot_reset as pr
import os

#pr.reset_plt(20,20)
#%%
tp_folder = '/scratch/s2555875/TP_files'

fig, ax = plt.subplots(tight_layout = True)
for file in os.listdir(tp_folder):
    tp = np.genfromtxt(os.path.join(tp_folder,file), skip_header = 1, comments = '#', names = True, dtype = None)
    ax.plot(tp['Temp'], tp['Pressure']/1e6, label = file[:6])
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('T [K]')
    ax.set_ylabel('P [bar]')
#fig.savefig('/scratch/s2555875/plot/helios_PTs_dist.pdf')
# %%
archean_Tp = np.genfromtxt('../atm/TP_helios.txt', dtype = None, skip_header=1, names=True, comments='#')
archean_Tp_0 = np.genfromtxt('/scratch/s2555875/HELIOS/output/2/2_tp.dat', dtype=None, skip_header=2, usecols=(1,2), names = ['Temp', 'Pressure'])
archean_Tp_1 = np.genfromtxt('/scratch/s2555875/HELIOS/output/3/3_tp.dat', dtype=None, skip_header=2, usecols=(1,2), names = ['Temp', 'Pressure'])
archean_Tp_2 = np.genfromtxt('/scratch/s2555875/HELIOS/output/archean/archean_tp.dat', dtype=None, skip_header=2, usecols=(1,2), names = ['Temp', 'Pressure'])
archean_Tp_3 = np.genfromtxt('/scratch/s2555875/HELIOS/output/archean_smooth/archean_smooth_tp.dat', dtype=None, skip_header=2, usecols=(1,2), names = ['Temp', 'Pressure'])
plt.axhline(0.14, color = 'k', linestyle = ':', label = 'tropopause')
plt.plot(archean_Tp['Temp'], archean_Tp['Pressure']*1e-6, label = 'X = 1e-42 above tropopause, used', linestyle = '-')
plt.plot(archean_Tp_0['Temp'], archean_Tp_0['Pressure']*1e-6, label = 'X = 1e-42 above tropopause', linestyle = '--')
plt.plot(archean_Tp_1['Temp'], archean_Tp_1['Pressure']*1e-6, label = 'X = 1e-99 above tropopause', linestyle = '-.')
plt.plot(archean_Tp_2['Temp'], archean_Tp_2['Pressure']*1e-6, label = 'X = 1e-6 above tropopause', linestyle = '-.')
#plt.plot(archean_Tp_3['Temp'], archean_Tp_3['Pressure']*1e-6, label = 'X = 1e-6 above tropopause, smooth', linestyle = '--')
plt.xlabel('T [K]')
plt.ylabel('P [bar]')
plt.yscale('log')
plt.gca().invert_yaxis()
plt.legend()
plt.savefig('/scratch/s2555875/plot/TP_comparison_waterprofile.pdf')
# %%
