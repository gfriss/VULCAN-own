import os
import subprocess # to run bash commands like in a terminal
import sys
import numpy as np
import parallel_functions as pf

nsim_dist = 15 # distances, for now
i_star = int(sys.argv[1]) # run number, used to keep track of which batch is being run (goes up to nsim_star)
scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
conv_file = os.path.join(scratch, 'converged.txt')
# ------setting up parameterspace for all runs------
# star type
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))
star_names = list(star_df.Name)
param_dict = {}
for star,a_min,a_max in zip(star_df.Name, star_df.a_min, star_df.a_max):
    distances = np.linspace(a_min, a_max, nsim_dist, endpoint = True)
    param_dict[star] = list(distances)
# defining network (crahcno is defualt), only used to identify sim folders and outputs
#network = ''
network = '_ncho'
# Boolian to check convergence and rerun if needed
check_conv = True
# ------end of parameter set up-----
for i_dist in range(nsim_dist): # since the star is given, need to only loop through distance cases
    sim = 'star_{}_'.format(star_names[i_star]) # param matrix first goes through the distances
    if i_dist < 10:
        sim += 'dist_0{}'.format(i_dist)                  # due to using this variable to allocate input files that are same for all the networks
        sim_folder = os.path.join(scratch, sim + network) # the network variable only comes into play when creating the sim folder and name
    else:
        sim += 'dist_{}'.format(i_dist)
        sim_folder = os.path.join(scratch, sim + network)
    # build files for simulation
    out_file = sim + network + '.vul'
    out_change = ','.join(['out_name', out_file, 'str'])
    output_dir_change = ','.join(['output_dir', output_folder, 'str'])
    new_cfg = os.path.join(sim_folder, 'vulcan_cfg.py')
    new_rad_file = pf.get_rad_prof(star_names[i_star])
    rad_file_change = ','.join(['sflux_file', new_rad_file, 'str'])
    new_r_star = str(star_df.loc[star_df.Name == star_names[i_star]].R.iloc[0])
    r_star_change = ','.join(['r_star', new_r_star, 'val'])
    new_orbit_radius = str(param_dict[star_names[i_star]][i_dist])
    orbit_radius_change = ','.join(['orbit_radius', str(new_orbit_radius), 'val'])
    new_tp_file = os.path.join(TP_folder, sim) + '.txt'
    tp_file_change = ','.join(['atm_file', new_tp_file, 'str'])
    # do a test on surface tmeperature so that only 0-100 surface tmeperature options will be simulated
    surface_temperatue = np.genfromtxt(new_tp_file, dtype = None, names = True, skip_header = 1, max_rows = 5)['Temp'][0]
    if surface_temperatue < 273 or surface_temperatue > 373:
        continue
    # first create simulation folder
    subprocess.check_call(['mkdir', sim_folder])
    # then make new cfg file
    subprocess.check_call(['python', 'gen_cfg.py', new_cfg, rad_file_change, r_star_change, orbit_radius_change, tp_file_change, out_change])
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
            out_file = sim + network + '_rerun.vul' # change this last so the initial composition will use the previous run
            out_change = ','.join(['out_name', out_file, 'str'])
            subprocess.check_call(['python', 'gen_cfg.py', new_cfg, out_change, vul_ini_change, ini_mix_change, yconv_min_change, 'rerun', sim])
            subprocess.check_call(['python', 'vulcan.py', '-n'])
    # then exit simulation folder and delete it
    os.chdir(wd)
    subprocess.check_call(['rm', '-rf', sim_folder])