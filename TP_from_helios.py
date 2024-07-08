import numpy as np
import os
import sys
import subprocess
from scipy import interpolate
#import parallel_functions as pf

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
    new_param_file = os.path.join(helios_folder, 'param_{}.dat'.format(sim_name))
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

    with open(os.path.join(scratch, 'TP_files', '{}_tp.txt'.format(sim_name)), 'w') as f:
        f.write('# (dyne/cm2) (K)     (cm2/s)\nPressure\tTemp\tKzz\n')
        for p,t,z in zip(P_vulcan,T_vulcan_helios,Kzz_vulcan):
            f.write('{:.3e}\t{:.1f}\t{:.3e}\n'.format(p,t,z))

#create_new_vulcan_pt('1')

a_list = np.linspace(0.723, 2, 15, endpoint = True) #HZ limits from Kopprapau et al. (2013) are 0.99 and 1.7, let's explore a bit more, from Venus to 2 au

for i,a in enumerate(a_list):
    sim = ''
    if i < 10:
        sim = 'sim_0{}_dist'.format(str(i))
    else:
        sim = 'sim_{}_dist'.format(str(i))

    new_p = create_param_file(a, sim)
    new_p = new_p
    #wd = os.getcwd()
    #os.chdir(helios_folder)
    subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p], cwd=helios_folder)
    #os.chdir(wd)
    create_new_vulcan_pt(sim)