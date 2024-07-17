#%%
import numpy as np
import os
import sys
import subprocess
from scipy import interpolate
#import parallel_functions as pf
#%%
scratch = '/scratch/s2555875'
helios_folder = os.path.join(scratch, 'HELIOS')
og_param = os.path.join(helios_folder, 'param.dat')
output = os.path.join(helios_folder, 'output')

#star_df = pf.read_stellar_data('/scratch/s2555875/stellar_flux/stellar_params.csv')
#star_names = star_df['Name'][:]

def create_param_file(dist, sim_name, sname = 'early_sun'):
    ''' Function to loop through the origianl HELIOS parameter file and change
        needed lines. For now these lines are: name and orbital distance.
        Future lines might contain: stellar spectrum file, stellar spectrum dataset,
        vertical mixing ratio'''
    new_str = ''
    with open(og_param, 'r') as og:
        for line in og.readlines():
            if line.strip() and line.split()[0] == 'name':
                bits = line.split()
                bits[2] = sim_name
                new_str += ' '.join(bits) + '\n'
            elif 'orbital distance' in line:
                bits = line.split()
                bits[6] = str(dist)
                new_str += ' '.join(bits) + '\n'
            else:
                new_str += line
    new_param_file = 'param_{}.dat'.format(sim_name)
    with open(new_param_file, 'w') as f:
        f.write(new_str)
    return new_param_file

def read_helios_tp(sim_name):
    tp_helios_file = os.path.join(helios_folder, 'output', sim_name, '{}_tp.dat'.format(sim_name))
    tp_helios_val = np.genfromtxt(tp_helios_file, dtype = None, skip_header=1, names = True, usecols = (1,2))
    T_helios = tp_helios_val['tempK']
    P_helios = tp_helios_val['press106bar']
    return T_helios, P_helios

def create_new_vulcan_pt(sim_name):
    T_h, P_h = read_helios_tp(sim_name)
    T_P_interp = interpolate.interp1d(P_h, T_h)

    vulcan_atm = np.genfromtxt('atm/T-P-Kzz_Pearce_B.txt', dtype = None, comments = '#', skip_header = 1, names = True)
    P_vulcan = vulcan_atm['Pressure']
    Kzz_vulcan = vulcan_atm['Kzz']
    T_vulcan_helios = T_P_interp(P_vulcan)

    with open(os.path.join(scratch, 'TP_files', '{}.txt'.format(sim_name)), 'w') as f:
        f.write('# (dyne/cm2) (K)     (cm2/s)\nPressure\tTemp\tKzz\n')
        for p,t,z in zip(P_vulcan,T_vulcan_helios,Kzz_vulcan):
            f.write('{:.3e}\t{:.1f}\t{:.3e}\n'.format(p,t,z))
#%%
create_new_vulcan_pt('3')
#%%
a_list = np.linspace(0.723, 2, 15, endpoint = True) #HZ limits from Kopprapau et al. (2013) are 0.99 and 1.7, let's explore a bit more, from Venus to 2 au

for i,a in enumerate(a_list):
    sim = ''
    if i < 10:
        sim = 'sim_0{}_dist'.format(str(i))
    else:
        sim = 'sim_{}_dist'.format(str(i))

    wd = os.getcwd()
    os.chdir(helios_folder)
    new_p = create_param_file(a, sim)
    subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p])#, cwd=helios_folder)
    os.chdir(wd)
    create_new_vulcan_pt(sim)
#%%
# creating fastchem styled mixing ratio file for HELIOS
location_early = '/scratch/s2555875/HELIOS/input/chemistry/early_earth/chem.dat'
location_current = '/scratch/s2555875/HELIOS/input/chemistry/earth/chem.dat'
mixing_file = 'atm/mixing_Pearce_B.txt'
def change_spec_name(sp_name):
    characters = []
    characters[:] = sp_name
    for i in range(len(characters)):
        if i < len(characters)-1:
            if characters[i+1].isnumeric() == False: # if another letter comes, it means there's only one of that element in the molecule
                characters[i] += '1'
        elif i == len(characters)-1:
            if characters[i].isnumeric() == False: # if another letter comes, it means there's only one of that element in the molecule
                characters[i] += '1'
    new_name = ''.join(characters)
    return new_name


def create_helios_vmr_from_vulcan(vul_vmr_file, vul_tp_file, filename):
    vul_vmr = np.genfromtxt(vul_vmr_file, dtype = None, comments = '#', skip_header = 1, names = True)
    vul_tp = np.genfromtxt(vul_tp_file, dtype = None, comments = '#', skip_header = 1, names = True)
    vmr_dict = {'P(bar)':vul_vmr['Pressure']/1e6, 'T(K)':vul_tp['Temp']}
    # keep ntot, ng and m as dummy now...
    vmr_dict['ntot'] = np.ones_like(vmr_dict['P(bar)'])
    vmr_dict['ng'] = np.ones_like(vmr_dict['P(bar)'])
    vmr_dict['m'] = np.ones_like(vmr_dict['P(bar)'])
    new_str = ''
    for name in vul_vmr.dtype.names:
        if name == 'Pressure':
            new_str += '#P(bar)\tT(k)\tn_<tot>(cm-3)\tn_g(cm-3)\tm(u)\t'
        else:
            new_str += '{}\t'.format(change_spec_name(name))
            vmr_dict[name] = vul_vmr[name]
    new_str += '\n'
    
    for i in range(len(vmr_dict['P(bar)'])):
        for key in vmr_dict.keys():
            new_str += '{:.4e}'.format(vmr_dict[key][i]) + '\t'
        new_str += '\n'

    with open(filename, 'w') as f:
        f.write(new_str)
    
# %%
