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
