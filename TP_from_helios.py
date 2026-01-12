#%%
import numpy as np
import os
import sys
import subprocess
from scipy import interpolate
import parallel_functions as pf
import pandas as pd

version = '_updated' # to use for updated initial mixing ratios, '' for old ones
scratch = '/scratch/s2555875'
helios_folder = os.path.join(scratch, 'HELIOS')
og_param = os.path.join(helios_folder, 'param{}.dat'.format(version))
output = os.path.join(helios_folder, 'output')
star_df = pd.read_csv(os.path.join(scratch, 'stellar_flux', 'stellar_params{}.csv'.format(version)))
nsim = 15
#%%
def change_line(line, bit_id, new_val):
    bits = line.split()
    bits[bit_id] = new_val
    return ' '.join(bits) + '\n'

def create_param_file(sim_name, dist = None, sname = None, manual = False, r_star = None, T_star = None, subfolder = ''):
    ''' Function to loop through the origianl HELIOS parameter file and change
        needed lines. For now these lines are: name and orbital distance.
        Future lines might contain: stellar spectrum file, stellar spectrum dataset,
        vertical mixing ratio'''
    new_str = ''
    with open(og_param, 'r') as og:
        for line in og.readlines():
            if line.strip() and line.split()[0] == 'name':
                new_str += change_line(line, 2, sim_name)
            elif 'realtime plotting' in line:
                new_str += change_line(line, 3, 'no')
            elif dist != None and 'orbital distance' in line:
                new_str += change_line(line, 6, str(dist))
            elif sname != None and 'path to stellar spectrum file' in line: # by default it asks for dile in param.dat
                new_str += change_line(line, 8, './star_tool/output/{}.h5'.format(sname.lower()))
            #elif sname != None and 'dataset in stellar spectrum file' in line: # by default it asks for dile in param.dat
            #    new_str += change_line(line, 8, '/r50_kdistr/phoenix/{}'.format(sname.lower()))
            elif sname != None and sname not in ['EARLY_SUN', 'SUN'] and 'dataset in stellar spectrum file' in line:
                new_str += change_line(line, 8, '/r50_kdistr/muscles/{}'.format(sname.lower()))
            elif sname != None and sname in ['EARLY_SUN', 'SUN'] and 'dataset in stellar spectrum file' in line:
                new_str += change_line(line, 8, '/r50_kdistr/ascii/{}'.format(sname.lower()))
            elif sname != None and 'planet =' in line:
                if manual:
                    new_str += change_line(line, 2, 'manual')
                else:
                    new_str += change_line(line, 2, sname.upper() + '_own') # all info on planets are in HELIOS/source/planet_database.py
            elif r_star != None and 'radius star [R_Sun]' in line:
                new_str += change_line(line, 6, str(r_star))
            elif T_star != None and 'temperature star [K]' in line:
                new_str += change_line(line, 6, str(T_star))
            else:
                new_str += line
    new_param_file = os.path.join(subfolder,'param_{}.dat'.format(sim_name))
    with open(new_param_file, 'w') as f:
        f.write(new_str)
    return new_param_file

def read_helios_tp(sim_name):
    tp_helios_file = os.path.join(helios_folder, 'output', sim_name, '{}_tp.dat'.format(sim_name))
    tp_helios_val = np.genfromtxt(tp_helios_file, dtype = None, skip_header=1, names = True, usecols = (1,2))
    T_helios = tp_helios_val['tempK']
    P_helios = tp_helios_val['press106bar']
    return T_helios, P_helios

def create_new_vulcan_pt(sim_name, subfolder = ''):
    T_h, P_h = read_helios_tp(sim_name)
    #T_P_interp = interpolate.interp1d(P_h, T_h, bounds_error = False, fill_value = "extrapolate")

    vulcan_atm = np.genfromtxt('/home/s2555875/VULCAN-2/atm/atm_Earth_Jan_Kzz.txt', dtype = None, comments = '#', skip_header = 1, names = True)
    P_vulcan = vulcan_atm['Pressure']
    Kzz_vulcan = vulcan_atm['Kzz']
    P_Kzz = interpolate.interp1d(P_vulcan, Kzz_vulcan, bounds_error = False, fill_value = "extrapolate")
    Kzz_h = P_Kzz(P_h)

    with open(os.path.join(scratch, 'TP_files', subfolder, '{}.txt'.format(sim_name)), 'w') as f:
        f.write('# (dyne/cm2) (K)     (cm2/s)\nPressure\tTemp\tKzz\n')
        for p,t,z in zip(P_h,T_h,Kzz_h):
            f.write('{:.3e}\t{:.1f}\t{:.3e}\n'.format(p,t,z))
#%%
#create_new_vulcan_pt('archean')
#%%
a_list = {'': np.linspace(0.839, 1.333, nsim, endpoint = True),
        '_updated': np.linspace(0.778, 1.307, nsim, endpoint = True)} # tested endpoints before running this cell to make sure durface temperature is habitable

def run_many_dist(dist_list):
    for i,a in enumerate(dist_list):
        sim = 'sim_{:02d}_dist{}'.format(i, version)

        wd = os.getcwd()
        os.chdir(helios_folder)
        new_p = create_param_file(sim, dist = a)
        subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p])#, cwd=helios_folder)
        os.chdir(wd)
        create_new_vulcan_pt(sim)

#run_many_dist(a_list[version])
#%%
def run_many_planets(star_table):
    name_list = star_table.Name[:]
    for i,name in enumerate(name_list):
        sim = 'sim_{:02d}_star{}'.format(i, version)
        
        wd = os.getcwd()
        os.chdir(helios_folder)
        new_p = create_param_file(sim, sname = name)
        subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p])#, cwd=helios_folder)
        os.chdir(wd)
        create_new_vulcan_pt(sim)
#run_many_planets(star_df)
#%%
def run_star_dist(star_table, factor = 1.1):
    param_matrix = []
    for star,a_min,a_max in zip(star_table.Name, star_table.a_min, star_table.a_max):
        #dist = pf.semi_major_list_from_Seff(star_df, star, nsim_dist, factor = factor)
        # for a few starts the original method got too cold surface temperatures so decided to test 
        # and make sure surface temperature are between 0and 100 Celsiusfor all
        dist = np.linspace(a_min, a_max, nsim, endpoint = True)
        for d in dist:
            param_matrix.append([star, d])

    for i,sim_i in enumerate(param_matrix):
        i_dist = i%nsim
        if i_dist in [0,14]:
            continue # extremities already done
        if sim_i[0] in ['EARLY_SUN']: # skipping sims that are done already
            subprocess.check_call(['cp', '-p', os.path.join(scratch, 'TP_files', 'sim_{:02d}_dist{}.txt'.format(i_dist, version)), os.path.join(scratch, 'TP_files', 'star_dist{}/star_{}_dist_{:02d}{}.txt'.format(version, sim_i[0], i_dist, version))])
            continue
        T = star_table.loc[star_table.Name == sim_i[0]].T_eff.iloc[0]
        R = star_table.loc[star_table.Name == sim_i[0]].R.iloc[0]
        sim = 'star_{}_dist_{:02d}{}'.format(sim_i[0], i_dist, version) # param matrix first goes through the distances
        wd = os.getcwd()
        os.chdir(helios_folder)
        new_p = create_param_file(sim, sname = sim_i[0], dist = sim_i[1], manual = True, T_star = T, r_star = R)
        subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p])#, cwd=helios_folder)
        os.chdir(wd)
        create_new_vulcan_pt(sim, subfolder = 'star_dist{}'.format(version))
run_star_dist(star_df)
# %%