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
# meteoritic bombardment
bomb_rate = np.linspace(3e23, 1e25, nsim) # values from Pearce et al. (2022) Fig A4 min and max
# for extra sims: np.linspace(3.5e23, 9e23, nsim)
vdep = '1.'#,1.e-4,0.,1.' # deposition velocity for H2, CO2, CH4, NH3 (CO from lighning too so have to adjust later...)
prod_sp = {'H2':0.65}#, 'CO2':1.32, 'CH4':1e-6, 'NH3':7e-5} # produced amount per impactor os m_mass from Zahnle et al (2020)
m_mass = 1e22 # stick to this for now by Zahnle et al. (2020)
# C/O ratio
co2_for_CtoO_range = np.linspace(0.001,0.1,nsim, endpoint = True)
# star type
star_df = pf.read_stellar_data(os.path.join(main_folder,'stellar_flux/stellar_params.csv')) # NOT necessarily nsim long...
a_star_list  = [0.0296, 0.0770, 0.0780, 0.1550, 0.1780, 0.3840, 0.4945, 0.6295, 0.6976, 1., 1.1623, 1.9047, 2.3155]
#phoenix_folder = os.path.join(main_folder, 'phoenix_fits')
#a_star_list = [0.0510, 0.0724, 0.0994, 0.1320, 0.1710, 0.3840, 0.4777, 0.5834, 0.7026, 0.8300, 0.9810, 1.1418, 2.5660, 2.9330, 3.3610]
#R_sol = 6.96e10 # cm
#Teff_list = np.array([2600+i*300 for i in range(nsim)])
#logg_list = np.array([5. for i in range(nsim)])
#logg_list[(Teff_list>=4000)&(Teff_list<6000)] = 4.5
#logg_list[Teff_list>=6000] = 4.0
#m = 0. # solar metallicity
# distance case
a_list = np.linspace(0.85, 1.35, nsim, endpoint = True) # tested endpoints before running this cell to make sure durface temperature is habitable
# in case new sims with HELIOS TP
helios_tp = ''
#helios_tp = 'helios_tp_'
# local meteorite case
h2_bar_list = np.linspace(0, 2, 15, endpoint = True)
# TOA pressure case
p_t_list = np.linspace(1e-2, 1e-1, 15, endpoint = True)
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
    sim = ''
    if i < 10:
        sim = helios_tp + 'sim_0' + str(i) + '_' + run_type # due to using this variable to allocate input files that are same for all the networks
        sim_folder = os.path.join(main_folder,sim + network) # the network variable only comes into play when creating the sim folder and name
    else:
        sim = helios_tp + 'sim_' + str(i) + '_' + run_type
        sim_folder = os.path.join(main_folder,sim + network)
    # build files for simulation
    out_file = sim + network + '.vul'
    out_change = 'out_name,' + out_file + ',str'
    new_cfg = os.path.join(sim_folder,'vulcan_cfg.py')
    # first create simulation folder
    subprocess.check_call(['mkdir', sim_folder])
    if run_type == 'meteor':
        # then BC because of meteorite and outgassing and lightning
        sp_names, sp_fluxes = pf.dict_to_input(pf.bombardment(bomb_rate[i], m_mass, prod_sp))
        # need: species, flux, vdep and name of txt file
        BC_bot_file = os.path.join(main_folder, 'BC_files', 'BC_bot_' + sim + '.txt')
        subprocess.check_call(['python', 'gen_BC.py', sp_names, sp_fluxes, vdep, BC_bot_file])
        # then use the new BC file in cfg file with output name
        BC_bot_change = ','.join(['bot_BC_flux_file', BC_bot_file, 'str'])
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, BC_bot_change, out_change])
    elif run_type == 'CtoO':
        # generate new mixing file
        new_mixing_file = os.path.join(main_folder, 'mixing_files', sim + 'mixing.txt')
        mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
        pf.gen_mixing(co2_for_CtoO_range[i], new_mixing_file)
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, out_change])
    elif run_type == 'star':
        # new stellar radiation files already created, need to update vulcan_cfg with filename and r_star and orbit_radius
        star_name = star_df.Name.iloc[i]
        #star_name = sim
        new_rad_file = pf.get_rad_prof(star_name)
        rad_file_change = ','.join(['sflux_file', new_rad_file, 'str'])
        new_r_star = str(star_df.loc[star_df.Name == star_name].R.iloc[0])
        #fits_file = os.path.join(phoenix_folder, '{:05d}_{:.2f}_{:.1f}.fits'.format(Teff_list[i],logg_list[i],m))
        #new_r_star = fits.open(fits_file)[0].header['PHXREFF'] / R_sol
        r_star_change = ','.join(['r_star', str(new_r_star), 'val'])
        new_orbit_radius = str(a_star_list[i])
        orbit_radius_change = ','.join(['orbit_radius', new_orbit_radius, 'val'])
        # new TP files are already created with HELIOS, need to update vulcan_cfg with filename and orbit_radius
        # surf temps are similar, but better to use corresponding ones
        tp_file = os.path.join(main_folder, 'TP_files', sim + '.txt')
        tp_change = ','.join(['atm_file', tp_file, 'str'])
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, rad_file_change, r_star_change, orbit_radius_change, tp_change, out_change])
    elif run_type == 'dist':
        # new TP files are already created with HELIOS, need to update vulcan_cfg with filename and orbit_radius
        tp_file = os.path.join(main_folder, 'TP_files', sim + '.txt')
        tp_change = ','.join(['atm_file', tp_file, 'str'])
        new_orbit_radius = str(a_list[i])
        orbit_radius_change = ','.join(['orbit_radius', new_orbit_radius, 'val'])
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, tp_change, orbit_radius_change, out_change])
    elif run_type == 'local':
        # generate new mixing ratios
        new_mixing_file = os.path.join(main_folder, 'mixing_files', sim + 'mixing.txt')
        mixing_change = ','.join(['vul_ini', new_mixing_file, 'str'])
        pf.gen_mixing_local(h2_bar_list[i], new_mixing_file)
        # H2 fiddles with cinvergence so reducing error tolerances
        atol_change = ','.join(['atol', str(1.E-7), 'val'])
        rtol_change = ','.join(['post_conden_rtol', str(0.5), 'val'])    
        # assumed to keep P-T profile and surface pressure the same
        # then change vulcan_cfg.py file
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, mixing_change, atol_change, rtol_change, out_change])
    elif run_type == 'pressure':
        # change top of atmosphere pressure
        p_t_change = ','.join(['P_t', str(p_t_list[i]), 'val'])
        subprocess.check_call(['python', 'gen_cfg.py', new_cfg, p_t_change, out_change])
    # then change to simulation folder and put symlinks in there to avoid copyying and make importing possible
    #subprocess.check_call(['cp', '-p', 'build_atm.py', 'chem_funs.py', 'op.py', 'phy_const.py', 'store.py', 'vulcan.py', sim_folder])
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
            out_file = sim + network + '_rerun.vul' # change this last so the initial composition will use the previous run
            out_change = ','.join(['out_name', out_file, 'str'])
            subprocess.check_call(['python', 'gen_cfg.py', new_cfg, out_change, vul_ini_change, ini_mix_change, yconv_min_change, 'rerun', sim])
            subprocess.check_call(['python', 'vulcan.py', '-n'])
    # then exit simulation folder and delete it
    os.chdir(wd)
    subprocess.check_call(['rm', '-rf', sim_folder])

# TO DO:
# DONE  need something to prepare atmosphere and find a good place to store that, chemistry will be the same
# DONE  similar for boundary conditions
# DONE  make a generalised, probably with inputs, script for initial mixing ratios
# DONE  write general vulcan_cfg generator -> NEED generic file to start from that has save structure
# DONE  do I need to copy vulcan.py for each function or can I call the same file? -> can call the same as they run as different processes with their own memory, variables, etc.
# DONE  after all these can the parallelisation run
# DONE  plus figure out the output structure: how to name and what groupings should I save them in
# DONE  so far every bit is in a different script, maybe can be as functions here, THINK about this -> less command line python stuff if they're here as functions...