import os
import subprocess # to run bash commands like in a terminal
import sys
import numpy as np
import parallel_functions as pf
from mpi4py.MPI import COMM_WORLD as CW # for paralellisation
import pickle

rank = CW.Get_rank()
size = CW.Get_size()

nsim = {'star': 13, 'dist': 15, 'CtoO': 15} # total number of simulations per type
nsim_total  = 1
for ns in nsim.values():
    nsim_total *= ns
sim_per_rank = int(nsim_total / size) + 1 # this is needed to distribure the tasks between the CPUs, +1 makes sure there will be enough CPUs to run all sims

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
conv_file = os.path.join(scratch, 'converged.txt')
non_conv_file = os.path.join(scratch, 'non_converged.txt')
# ------setting up parameterspace for all runs------
# star type
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))
star_names = list(star_df.Name)
T_eff = list(star_df.T_eff)

C_to_O = []
for j in range(nsim['CtoO']):
    sim_CtoO = 'sim_'
    if j < 10:
        sim_CtoO += '0{}_CtoO'.format(j)
    else:
        sim_CtoO += '{}_CtoO'.format(j)
    mixing_file = os.path.join(scratch, 'mixing_files', sim_CtoO + 'mixing.txt')
    with open(os.path.join(scratch, 'output', sim_CtoO + '.vul'), 'rb') as handle:
        data_CtoO = pickle.load(handle)
    C_to_O.append(pf.calc_C_to_O(data_CtoO, mixing_file))        
        
# parameter dictionary of stars and distances, C/O is not inlcuded because it's same all over
param_dict = {}
for star,a_min,a_max in zip(star_df.Name, star_df.a_min, star_df.a_max):
    distances = np.linspace(a_min, a_max, nsim['dist'], endpoint = True)
    param_dict[star] = list(distances)
# defining network (crahcno is defualt), only used to identify sim folders and outputs
#network = ''
network = '_ncho'
# Boolian to check convergence and rerun if needed
check_conv = True
# ------end of parameter set up-----
# ------start of simulation loop------
for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):   # paralellisation itself, it spreads the task between the CPUs
                                                            # this is the magic, after this just think of it as a normal, sequential loop
    if i > nsim_total - 1: # in case using total sim number is not divisible by number of CPUs
        continue
    i_star = i//(nsim['dist']*nsim['CtoO']) # which star (changes after looped through all distance and C/O possibilities)
    i_dist = (i//nsim['CtoO'])%nsim['dist'] # which distance (changes after looped through all C/O possibilities and then restarts when all distances are done)
    i_CtoO = i%nsim['CtoO'] # which C/O (simply loops through the C/O possibilities)
    # build simulation and folder names
    sim = 'star_{}'.format(star_names[i_star]) # param matrix first goes through the stars
    dist_sim = 'dist_' + pf.get_str_number(i_dist) # then the distance
    CtoO_sim = 'CtoO_' + pf.get_str_number(i_CtoO) # then the C/O ratio
    TP_sim = '_'.join([sim, dist_sim]) # TP profile is the same for all C/O ratios
    sim = '_'.join([sim, dist_sim, CtoO_sim]) # this is the final simulation name
    sim_folder = os.path.join(scratch, sim + network)
    # skip parts that are already done, i.e. all distance with base C/O (0th) (same name as used for TP profile)
    # but add their output to the file
    # to be deleted later
    #with open(conv_file, 'r') as f:
    #    conv_text = f.read()
    #if sim + network + '.vul' in conv_text:
    #    continue
    #with open(non_conv_file, 'r') as f:
    #    non_conv_text = f.read()
    #if sim+network+'_rerun.vul' in non_conv_text:
    #    continue # don't do it again if it didn't converge the first time
    #if i_CtoO == 0:
    #    with open(os.path.join(output_folder, TP_sim+network+'.vul'), 'rb') as handle:
    #        data = pickle.load(handle)
    #    hcn_rain = pf.rainout(data, rain_spec = 'HCN_rain', g_per_mol = 27)
    #    h20_rain = pf.rainout(data, rain_spec = 'H2O_rain', g_per_mol = 18)
    #    with open(os.path.join(scratch, 'star_dist_CtoO_rain.txt'), 'a') as f:
    #        f.write('{}\t{}\t{}\t{}\t{}\n'.format(T_eff[i_star], param_dict[star_names[i_star]][i_dist], C_to_O[i_CtoO], hcn_rain, h20_rain))
    #    continue
    # build files for simulation
    out_file = sim + network + '_nowash.vul'
    out_change = ','.join(['out_name', out_file, 'str'])
    output_dir_change = ','.join(['output_dir', output_folder, 'str'])
    new_cfg = os.path.join(sim_folder, 'vulcan_cfg.py')
    new_rad_file = pf.get_rad_prof(star_names[i_star])
    rad_file_change = ','.join(['sflux_file', new_rad_file, 'str'])
    new_r_star = str(star_df.loc[star_df.Name == star_names[i_star]].R.iloc[0])
    r_star_change = ','.join(['r_star', new_r_star, 'val'])
    new_orbit_radius = str(param_dict[star_names[i_star]][i_dist])
    orbit_radius_change = ','.join(['orbit_radius', str(new_orbit_radius), 'val'])
    new_tp_file = os.path.join(TP_folder, TP_sim+'.txt')
    tp_file_change = ','.join(['atm_file', new_tp_file, 'str'])
    new_mixing_file = os.path.join(scratch, 'mixing_files', 'sim_{}_CtoOmixing.txt'.format(pf.get_str_number(i_CtoO)))
    mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
    # first create simulation folder
    subprocess.check_call(['mkdir', sim_folder])
    # then make new cfg file
    subprocess.check_call(['python', 'gen_cfg.py', new_cfg, rad_file_change, r_star_change, orbit_radius_change, tp_file_change, mixing_change, out_change])
    # then change to simulation folder and put symlinks in there to avoid copyying and make importing possible
    wd = os.getcwd()
    os.chdir(sim_folder)
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/build_atm.py', 'build_atm.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/make_chem_funs.py', 'make_chem_funs.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/op.py', 'op.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/phy_const.py', 'phy_const.py'])
    subprocess.check_call(['cp', '-p', '/home/s2555875/VULCAN-2/store.py', 'store.py'])
    subprocess.check_call(['cp', '-p', '/home/s2555875/VULCAN-2/vulcan.py', 'vulcan.py'])
    subprocess.check_call(['cp', '-r', '/home/s2555875/VULCAN-2/thermo', 'thermo'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/atm', 'atm'])
    subprocess.check_call(['cp', '-p', '/home/s2555875/VULCAN-2/gen_cfg.py', 'gen_cfg.py'])
    # then run vulcan.py
    subprocess.check_call(['python', 'vulcan.py'])
    # then check convergence and rerun once if needed:
    if check_conv:
        with open(conv_file, 'r') as f:
            conv_text = f.read()
        if out_file not in conv_text:
            vul_ini_change = ','.join(['vul_ini', os.path.join(output_folder,out_file), 'str'])
            ini_mix_change = ','.join(['ini_mix', 'vulcan_ini', 'str'])
            yconv_min_change = ','.join(['yconv_min', str(0.2), 'val']) # 0.1 is the default, allow double
            out_file = sim + network + '_nowash_rerun.vul' # change this last so the initial composition will use the previous run
            out_change = ','.join(['out_name', out_file, 'str'])
            subprocess.check_call(['python', 'gen_cfg.py', new_cfg, out_change, vul_ini_change, ini_mix_change, yconv_min_change, 'rerun', sim])
            subprocess.check_call(['python', 'vulcan.py', '-n'])
    # save results
    hcn_rain, h20_rain = np.nan, np.nan # will keep NaN if not converged
    with open(conv_file, 'r') as f:
        conv_text = f.read()
    if out_file in conv_text:
        with open(os.path.join(output_folder, out_file), 'rb') as handle:
            data = pickle.load(handle)
        hcn_rain = pf.rainout(data, rain_spec = 'HCN_rain', g_per_mol = 27)
        h20_rain = pf.rainout(data, rain_spec = 'H2O_rain', g_per_mol = 18)
    with open(os.path.join(scratch, 'star_dist_CtoO_rain.txt'), 'a') as f:
        f.write('{}\t{}\t{}\t{}\t{}\n'.format(T_eff[i_star], param_dict[star_names[i_star]][i_dist], C_to_O[i_CtoO], hcn_rain, h20_rain))
    # then exit simulation folder and delete it along with output and cfg file
    os.chdir(wd)
    subprocess.check_call(['rm', '-rf', sim_folder])
    subprocess.check_call(['rm', os.path.join(output_folder, out_file)])
    subprocess.check_call(['rm', os.path.join(output_folder, 'cfg_'+out_file[:-4]+'.txt')])