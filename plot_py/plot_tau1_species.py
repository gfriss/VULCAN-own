import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle
       
nsim = 15
version = ''
plot_name_list = ['sim_{:02d}_CH4-balanced{}_tau'.format(i, version) for i in range(nsim)] # ['archean_ncho_nowash{}_tau.format(version)] or ['sim_{:02d}_dist{}_tau'.format(i, version) for i in range(nsim)]
vul_data_list = ['/scratch/s2555875/output/sim_{:02d}_CH4-balanced_ncho_nowash{}.vul'.format(i, version) for i in range(nsim)] # ['/scratch/s2555875/output/archean{}_ncho_nowash.vul'.format(version)] or ['/scratch/s2555875/output/sim_{:02d}_dist_ncho_nowash{}.vul'.format(i, version) for i in range(nsim)]
plot_dir = '/scratch/s2555875/plot/tau/'

colors = ['c','b','g','r','m','y','chocolate','orange','pink','grey','darkred','salmon','steelblue','hotpink','k','c','b','g','r','m','y']

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$','O3':'O$_3$' }
photo_sp = ['N2', 'CO2', 'CO', 'CH4', 'H2O', 'C2H4', 'C2H2','C2H6']
scat_sp = ['N2', 'CO2']
# photosphere
tau1 = 1.
for plot_name, vul_data in zip(plot_name_list, vul_data_list):
    color_index = 0
    vul_rerun = vul_data.replace('.vul', '_rerun.vul')
    if os.path.exists(vul_rerun):
        vul_data = vul_rerun
    with open(vul_data, 'rb') as handle:
        data = pickle.load(handle)
        
    # photosphere of each species
    tau_sp= {}

    bins = data['variable']['bins']
    dz = data['atm']['dz']
    nz = len(dz)
    vulcan_spec = data['variable']['species']

    tau_sum = np.zeros( (nz+1, len(bins)) )

    tau_sp, tau_scat = {}, {}

    for sp in photo_sp: tau_sp[sp] = np.zeros( (nz+1, len(bins)) )
    tau_scat = np.zeros( (nz+1, len(bins)) )
                
    for j in range(nz-1,-1,-1):
        for sp in photo_sp: # scat_sp are not necessary photo_sp, e.g. He
            #absorption of speecies sp at level j
            tau_sp[sp][j] += data['variable']['y'][j,vulcan_spec.index(sp)] * dz[j] * data['variable']['cross'][sp]
            #absorption of all speecies sp at level j
            tau_sum[j] += tau_sp[sp][j]
            
            #adding the layer above for species j
            tau_sp[sp][j] += tau_sp[sp][j+1]
        
        for sp in scat_sp:
            tau_scat[j] += data['variable']['y'][j,vulcan_spec.index(sp)] * dz[j] * data['variable']['cross_scat'][sp]
            tau_sum[j] += data['variable']['y'][j,vulcan_spec.index(sp)] * dz[j] * data['variable']['cross_scat'][sp]
        
        # adding the layer above only at the end of species loop    
        tau_scat[j] += tau_scat[j+1] 
        tau_sum[j] += tau_sum[j+1] 

        
    plt.figure()


    photosph = []
    photosph_abs = []
    photosph_sp = {}
    photosph_scat = []

    photosph_sum = []

    for n in range(len(bins)):
        photosph.append( data['atm']['pico'][np.argmin(np.abs(data['variable']['tau'][:,n]-tau1)) ]/1.e5 )
        #photosph_abs.append( data['atm']['pco'][np.argmin(np.abs(data['variable']['tau_abs'][:,n]-tau1)) ]/1.e6 )
        photosph_scat.append( data['atm']['pico'][np.argmin(np.abs(tau_scat[:,n]-tau1)) ]/1.e5 )
        
        photosph_sum.append( data['atm']['pico'][np.argmin(np.abs(tau_sum[:,n]-tau1)) ]/1.e5 )
        
    for sp in photo_sp:
        photosph_sp[sp] = []
        for n in range(len(bins)):
            photosph_sp[sp].append( data['atm']['pico'][np.argmin(  np.abs(tau_sp[sp][:,n] -tau1)) ]/1.e5 )
        
        



    plt.plot(bins, photosph, c='k', lw=1.5, alpha=0.7, label=r'$\tau$=1 (total)')
    for sp in photo_sp:
        # if sp == 'H2O':
    #         plt.plot(bins, photosph_sp[sp], color=colors[color_index], label = tex_labels[sp], alpha=0.6, lw=2 )
    #     else:
        if sp in tex_labels: tex_lab = tex_labels[sp]
        else: tex_lab = sp
        plt.plot(bins, photosph_sp[sp], color=colors[color_index], label = tex_lab, alpha=0.6, lw=2 )
        color_index += 1

    plt.plot(bins, photosph_scat, color=colors[color_index+1], label = 'Rayleigh', alpha=0.5 , lw=1.2)


    #plt.plot(bins, photosph_sum, c='plum', lw=3, alpha=0.6)


    #plt.plot(data['variable']['bins'], photosph_abs, c='red')
            
    plt.gca().set_yscale('log') 
    plt.gca().invert_yaxis() 
    #plt.ylim((0,95))
    plt.xlim((0,300))
    plt.legend(frameon=0, prop={'size':13}, loc='best')
    plt.xlabel("wavelength (nm)", fontsize=16)
    plt.ylabel("Pressure (bar)", fontsize=16)
    plt.title('Photosphere - ' + plot_name[:-4])
    plt.savefig(plot_dir + plot_name + '.pdf')
    plt.close()