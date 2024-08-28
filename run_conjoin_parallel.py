import os
import subprocess # to run bash commands like in a terminal
import sys
import numpy as np
import parallel_functions as pf
from mpi4py.MPI import COMM_WORLD as CW # for paralellisation

rank = CW.Get_rank()
size = CW.Get_size()

nsim_dist = 15 # distances, for now
nsim_star = 13 # stars, for now
sim_per_rank = int(nsim_dist*nsim_star / size) + 1 # this is needed to distribure the tasks between tha CPUs

scratch = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(scratch, 'output/star_dist/')
TP_folder = os.path.join(scratch, 'TP_files/star_dist')
conv_file = os.path.join(scratch, 'converged.txt')
check_conv = True
# ------setting up parameterspace for all runs------
# star type
star_df = pf.read_stellar_data(os.path.join(scratch, 'stellar_flux/stellar_params.csv'))
param_matrix = []
for star in star_df.Name:
    dist = pf.semi_major_list_from_Seff(star_df, star, nsim_dist, factor = 1.1)
    for d in dist:
        param_matrix.append([star, d])

# ------end of parameter set up-----
for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):   # paralellisation itself, it spreads the task between the CPUs
                                                            # this is the magic, after this just think of it as a normal, sequential loop
    if i > nsim_dist*nsim_star: # in case using total sim number is not divisible by number of CPUs
        continue
    i_star = i//nsim_dist
    i_dist = i%nsim_dist
    sim = 'star_{}_'.format(param_matrix[i][0]) # param matrix first goes through the distances
    if i_dist < 10:
        sim += 'dist_0{}'.format(i_dist)
        sim_folder = os.path.join(scratch, sim)
    else:
        sim += 'dist_{}'.format(i_dist)
        sim_folder = os.path.join(scratch, sim)
    # build files for simulation
    out_file = sim + '.vul'
    out_change = ','.join(['out_name', out_file, 'str'])
    output_dir_change = ','.join(['output_dir', output_folder, 'str'])
    new_cfg = os.path.join(sim_folder, 'vulcan_cfg.py')
    new_rad_file = pf.get_rad_prof(param_matrix[i][0])
    rad_file_change = ','.join(['sflux_file', new_rad_file, 'str'])
    new_r_star = str(star_df.loc[star_df.Name == param_matrix[i][0]].R.iloc[0])
    r_star_change = ','.join(['r_star', new_r_star, 'val'])
    new_orbit_radius = str(param_matrix[i][1])
    orbit_radius_change = ','.join(['orbit_radius', str(new_orbit_radius), 'val'])
    new_tp_file = os.path.join(TP_folder, sim) + '.txt'
    tp_file_change = ','.join(['atm_file', new_tp_file, 'str'])
    # do a test on surface tmeperature so that only 0-100 surface tmeperature options will be simulated
    surface_temperatue = np.genfromtxt(new_tp_file, dtype = None, names = True, skip_header = 1, max_rows = 5)['Temp'][0]
    if surface_temperatue < 273 or surface_temperatue > 373:
        continue
    # do test on convergence and rerun what didn't converge wiht a changed atol and rtol
    if check_conv:
        with open(conv_file, 'r') as f:
            conv_text = f.read()
        if out_file not in conv_text:
            atol_change = ','.join(['atol', str(1.E-3), 'val'])
            rtol_change = ','.join(['post_conden_rtol', str(0.1), 'val'])
            # first create simulation folder
            subprocess.check_call(['mkdir', sim_folder])
            # then make new cfg file
            subprocess.check_call(['python', 'gen_cfg.py', new_cfg, rad_file_change, r_star_change, orbit_radius_change, tp_file_change, out_change, atol_change, rtol_change])
        else:
            continue
    else:
        # first create simulation folder
        subprocess.check_call(['mkdir', sim_folder])
        # then make new cfg file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, rad_file_change, r_star_change, orbit_radius_change, tp_file_change, out_change])
    # then change to simulation folder and put symlinks in there to avoid copyying and make importing possible
    wd = os.getcwd()
    os.chdir(sim_folder)
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/build_atm.py', 'build_atm.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/chem_funs.py', 'chem_funs.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/op.py', 'op.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/phy_const.py', 'phy_const.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/store.py', 'store.py'])
    subprocess.check_call(['cp', '-p', '/home/s2555875/VULCAN-2/vulcan.py', 'vulcan.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/thermo', 'thermo'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/atm', 'atm'])
    
    # then run vulcan.py
    subprocess.check_call(['python', 'vulcan.py', '-n'])
    # then exit simulation folder and delete it
    os.chdir(wd)
    subprocess.check_call(['rm', '-rf', sim_folder])