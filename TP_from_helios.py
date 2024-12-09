#%%
import numpy as np
import os
import sys
import subprocess
from scipy import interpolate
import parallel_functions as pf
import pandas as pd

scratch = '/scratch/s2555875'
helios_folder = os.path.join(scratch, 'HELIOS')
og_param = os.path.join(helios_folder, 'param.dat')
output = os.path.join(helios_folder, 'output')
star_df = pd.read_csv(os.path.join(scratch, 'stellar_flux', 'stellar_params.csv'))
nsim_dist = 15
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
a_list = np.linspace(0.85, 1.35, nsim_dist, endpoint = True) # tested endpoints before running this cell to make sure durface temperature is habitable

def run_many_dist(dist_list):
    for i,a in enumerate(dist_list):
        sim = ''
        if i < 10:
            sim = 'sim_0{}_dist'.format(str(i))
        else:
            sim = 'sim_{}_dist'.format(str(i))

        wd = os.getcwd()
        os.chdir(helios_folder)
        new_p = create_param_file(sim, dist = a)
        subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p])#, cwd=helios_folder)
        os.chdir(wd)
        create_new_vulcan_pt(sim)

#run_many_dist(a_list)
#%%
def run_many_planets(star_table):
    name_list = star_table.Name[:]
    for i,name in enumerate(name_list):
        sim = ''
        if i < 10:
            sim = 'sim_0{}_star'.format(str(i))
        else:
            sim = 'sim_{}_star'.format(str(i))
        
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
        dist = np.linspace(a_min, a_max, nsim_dist, endpoint = True)
        for d in dist:
            param_matrix.append([star, d])

    for i,sim_i in enumerate(param_matrix):
        i_dist = i%nsim_dist
        #if i_dist not in [0,14]:
        #    continue # testing extremities
        #if sim_i[0] in ['EARLY_SUN', 'SUN', 'GJ1214', 'GJ674', 'GJ15A', 'HD40307', 'TRAPPIST-1', 'GJ676A', 'HD97658', 'HD149026', 'WASP17','GJ729']:
        #    continue
        #if sim_i[0] in ['HD85512'] and i_dist == 14:
        #    continue
        T = star_table.loc[star_table.Name == sim_i[0]].T_eff.iloc[0]
        R = star_table.loc[star_table.Name == sim_i[0]].R.iloc[0]
        sim = 'star_{}_'.format(sim_i[0]) # param matrix first goes through the distances
        if i_dist < 10:
            sim += 'dist_0{}'.format(i_dist)
        else:
            sim += 'dist_{}'.format(i_dist)
        wd = os.getcwd()
        os.chdir(helios_folder)
        new_p = create_param_file(sim, sname = sim_i[0], dist = sim_i[1], manual = True, T_star = T, r_star = R)
        subprocess.check_call(['python', 'helios.py', '-parameter_file', new_p])#, cwd=helios_folder)
        os.chdir(wd)
        create_new_vulcan_pt(sim, subfolder = 'star_dist')
run_star_dist(star_df)
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


def create_fastchem_like_vmr_from_vulcan(vul_vmr_file, vul_tp_file, filename):
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
import pickle
def create_helios_vmr_from_vulcan(vul_file, filename, species):
    with open(vul_file, 'rb') as handle:
        data = pickle.load(handle)
    spec_list = data['variable']['species']
    pressure = data['atm']['pco'] # in cgs
    new_str = '# (dyne/cm2)\nPressure\t'
    vmrs = {'Pressure': pressure}
    for sp in species:
        vmrs[sp] = data['variable']['ymix'][:, spec_list.index(sp)]
        new_str += '{}\t'.format(sp)
    
    new_str += '\n'
    for i in range(len(pressure)):
        for k in vmrs.keys():
            new_str += '{:.3e}\t'.format(vmrs[k][i])
        new_str += '\n'

    with open(filename, 'w') as f:
        f.write(new_str)

# %%
