import os
import subprocess # to run bash commands like in a terminal
import sys
import numpy as np
import parallel_functions as pf
from astropy.io import fits
from mpi4py.MPI import COMM_WORLD as CW # for paralellisation

run_type = sys.argv[1] # 'meteor' or 'CtoO' or 'star' according to which experiment we run, used to make loop at end more readable

rank = CW.Get_rank()
size = CW.Get_size()

nsim = 15
sim_per_rank = int(nsim / size) # this is needed to distribure the tasks between tha CPUs

main_folder = '/scratch/s2555875' # place to store outputs
output_folder = os.path.join(main_folder,'output')
conv_file = os.path.join(main_folder, 'converged.txt')
check_conv = True
# ------setting up parameterspace for all runs------
# version (nowash and for updated methane conc.)
nowash = '_nowash'
version = '_updated' # '_updated' or '' for old
mixing_file = {'': './atm/mixing_table_archean.txt', '_updated': './atm/mixing_table_archean_updated.txt'}
# meteoritic bombardment
bomb_rate = np.linspace(3e23, 1e25, nsim) # values from Pearce et al. (2022) Fig A4 min and max
# for extra sims: np.linspace(3.5e23, 9e23, nsim)
vdep = '1.'#,1.e-4,0.,1.' # deposition velocity for H2, CO2, CH4, NH3 (CO from lighning too so have to adjust later...)
prod_sp = {'H2':0.65}#, 'CO2':1.32, 'CH4':1e-6, 'NH3':7e-5} # produced amount per impactor os m_mass from Zahnle et al (2020)
m_mass = 1e22 # stick to this for now by Zahnle et al. (2020)
# C/O ratio
co2_for_CtoO_range = np.linspace(0.1, 0.001, nsim, endpoint = True)
# star type
star_df = pf.read_stellar_data(os.path.join(main_folder,'stellar_flux/stellar_params{}.csv'.format(version)))
a_star_list  = {'': [0.0269, 0.0697, 0.0719, 0.1451, 0.1679, 0.3653, 0.4808, 0.6228, 0.6960, 1., 1.1643, 1.8887, 2.3367], # NOT necessarily nsim long...
                '_updated': [0.0296, 0.0757, 0.0761, 0.1524, 0.1762, 0.3780, 0.4915, 0.6304, 0.7027, 1., 1.1713, 1.8852, 2.3277]}
# distance case
a_list = {'': np.linspace(0.839, 1.333, nsim, endpoint = True),
        '_updated': np.linspace(0.778, 1.307, nsim, endpoint = True)} # tested endpoints before running this cell to make sure surface temperature is habitable
# methane test case
ch4_range = np.logspace(-6, np.log10(5e-3), 15, endpoint = False) # mixing ratio range for methane, 1-5000ppmv (ommit endpoint to avoid overlap with Achean case)
# defining network (crahcno is defualt), only used to identify sim folders and outputs
#network = ''
network = '_ncho'
# Boolian to check convergence and rerun if needed
check_conv = True

# ------end of parameter set up-----
for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):   # paralellisation itself, it spreads the task between the CPUs
                                                            # this is the magic, after this just think of it as a normal, sequential loop
    if run_type == 'star' and i >= len(star_df.Name): # fewer stars than 15...
        continue
    sim = 'sim_{:02d}_{}'.format(i, run_type) # input files use this and they are the same for all versions and networks
    sim_folder = os.path.join(main_folder,sim + network + version)
    # build files for simulation
    out_file = sim + network + nowash + version + '.vul'
    out_change = 'out_name,' + out_file + ',str'
    new_cfg = os.path.join(sim_folder,'vulcan_cfg.py')
    mixing_change = ','.join(['vul_ini', mixing_file[version], 'str'])
    tp_file = os.path.join(main_folder, 'TP_files', '/scratch/s2555875/TP_files/archean{}.txt'.format(version))
    tp_change = ','.join(['atm_file', tp_file, 'str'])
    # first create simulation folder
    subprocess.check_call(['mkdir', sim_folder])
    if run_type == 'meteor':
        # then BC because of meteorite and outgassing and lightning
        sp_names, sp_fluxes = pf.dict_to_input(pf.bombardment(bomb_rate[i], m_mass, prod_sp))
        BC_bot_file = os.path.join(main_folder, 'BC_files', 'BC_bot_' + sim + '.txt')
        subprocess.check_call(['python', 'gen_BC.py', sp_names, sp_fluxes, vdep, BC_bot_file])
        # then use the new BC file in cfg file with output name
        BC_bot_change = ','.join(['bot_BC_flux_file', BC_bot_file, 'str'])
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, tp_change, BC_bot_change, out_change])
    elif run_type == 'CtoO':
        # generate new mixing file
        new_mixing_file = os.path.join(main_folder, 'mixing_files', sim + version + 'mixing.txt')
        mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
        pf.gen_mixing(co2_for_CtoO_range[i], new_mixing_file, version = version)
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, tp_change, out_change])
    elif run_type == 'star':
        # new stellar radiation files already created, need to update vulcan_cfg with filename and r_star and orbit_radius
        star_name = star_df.Name.iloc[i]
        new_rad_file = pf.get_rad_prof(star_name)
        rad_file_change = ','.join(['sflux_file', new_rad_file, 'str'])
        new_r_star = str(star_df.loc[star_df.Name == star_name].R.iloc[0])
        r_star_change = ','.join(['r_star', str(new_r_star), 'val'])
        new_orbit_radius = str(a_star_list[version][i])
        orbit_radius_change = ','.join(['orbit_radius', new_orbit_radius, 'val'])
        # new TP files are already created with HELIOS, need to update vulcan_cfg with filename and orbit_radius
        # surf temps are similar, but better to use corresponding ones
        tp_file = os.path.join(main_folder, 'TP_files', sim + version + '.txt')
        tp_change = ','.join(['atm_file', tp_file, 'str'])
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, rad_file_change, r_star_change, orbit_radius_change, tp_change, out_change])
    elif run_type == 'dist':
        # new TP files are already created with HELIOS, need to update vulcan_cfg with filename and orbit_radius
        tp_file = os.path.join(main_folder, 'TP_files', sim + version + '.txt')
        tp_change = ','.join(['atm_file', tp_file, 'str'])
        new_orbit_radius = str(a_list[version][i])
        orbit_radius_change = ','.join(['orbit_radius', new_orbit_radius, 'val'])
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, tp_change, orbit_radius_change, out_change])
    elif run_type == 'methane':
        # generate new mixing file
        new_mixing_file = os.path.join(main_folder, 'mixing_files', sim + 'mixing.txt')
        mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
        pf.gen_mixing(ch4_range[i], new_mixing_file, species = 'CH4') # convert bar to dyne/cm2
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, out_change])
    elif run_type == 'CH4-H2':
        # generate new mixing file
        new_mixing_file = os.path.join(main_folder, 'mixing_files', sim + 'mixing.txt')
        mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
        pf.gen_mixing(ch4_range[i], new_mixing_file, species = 'CH4-H2') # convert bar to dyne/cm2
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, out_change])
    elif run_type == 'CH4-balanced':
        # generate new mixing file
        new_mixing_file = os.path.join(main_folder, 'mixing_files', sim + 'mixing.txt')
        mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
        pf.gen_mixing(ch4_range[i], new_mixing_file, species = 'CH4-balanced') # convert bar to dyne/cm2
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, out_change])
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
    # check convergence and rerun once if needed:
    if check_conv:# == True and n_round == 1:
        with open(conv_file, 'r') as f:
            conv_text = f.read()
        if out_file not in conv_text:
            vul_ini_change = ','.join(['vul_ini', os.path.join(output_folder,out_file), 'str'])
            ini_mix_change = ','.join(['ini_mix', 'vulcan_ini', 'str'])
            yconv_min_change = ','.join(['yconv_min', str(0.2), 'val']) # 0.1 is the default, allow double
            out_file = sim + network + nowash + version +  '_rerun.vul' # change this last so the initial composition will use the previous run
            out_change = ','.join(['out_name', out_file, 'str'])
            subprocess.check_call(['python', 'gen_cfg.py', new_cfg, out_change, vul_ini_change, ini_mix_change, yconv_min_change, 'rerun', sim])
            subprocess.check_call(['python', 'vulcan.py', '-n'])
    # then exit simulation folder and delete it
    os.chdir(wd)
    subprocess.check_call(['rm', '-rf', sim_folder])