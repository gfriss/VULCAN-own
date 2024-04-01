import os
import subprocess # to run bash commands like in a terminal
import numpy as np
#from phy_const import Navo, kb
# this is how to use it:
# subprocess.run(['everything', 'you', 'would', 'type', 'word', 'by', 'word'], check = True, text = True)
from mpi4py.MPI import COMM_WORLD as CW # for paralellisation

#from phy_const import Navo, kb    
Navo = 6.02214086E+23
# BC from meteorite bombardment rate:
def bar_to_number(x): # 100 bar is 2300 mol/cm-2 according to Zahnler et al. (2020)
    return x * 2300 * Navo / 100

def gpergyr_to_gpers(x):
    return x / 3.154e16

def bombardment(rate, mass, prod):
    ''' Function to calculate the surface emission due to meteoritic bombardment using Zahnler et al. (2020)
        equilibrium production (Table1 last rows). Rate is in g/Gyr, mass is in g and prod is a dict
        containing the molecules and their eq. production in bars.'''
    bc = {}
    for sp in prod.keys():
        nd_per_mass = bar_to_number(prod[sp]) / mass # 1/(cm2 g)
        surface_flux = nd_per_mass * gpergyr_to_gpers(rate)
        bc[sp] = str(surface_flux)
    return bc

def dict_to_input(d):
    k = ','.join(d.keys())
    v = []
    for val in d.values():
        v.append(val)
    val_str = ','.join(v)
    return k, val_str

rank = CW.Get_rank()
size = CW.Get_size()

nsim = 15
sim_per_rank = int(nsim / size) # this is needed to distribure the tasks between tha CPUs

main_folder = '/scratch/s2555875/' # place to store outputs
bomb_rate = np.linspace(3e23, 1e25, nsim) # values from Pearce et al. (2022) Fig A4 min and max
vdep = '1.,1.e-4,0.,1.' # deposition velocity for H2, CO2, CH4, NH3 (CO from lighning too so have to adjust later...)
prod_sp = {'H2':0.65, 'CO2':1.32, 'CH4':1e-6, 'NH3':7e-5} # produced amount per impactor os m_mass from Zahnle et al (2020)
m_mass = 1e22 # stick to this for now by Zahnle et al. (2020)

for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):   # paralellisation itself, it spreads the task between the CPUs
                                                            # this is the magic, after this just think of it as a normal, sequential loop
    sim = ''
    if i < 10:
        sim = 'sim_0' + str(i)
        sim_folder = main_folder + sim
        #new_folder = os.path.join(main_folder, sim) # my folder naming conventions
    else:
        sim = 'sim_' + str(i)
        sim_folder = main_folder + sim
        #new_folder = os.path.join(main_folder, sim)
    # build files for simulation
    # first create simulation folder
    subprocess.check_call(['mkdir', sim_folder])
    # then BC because of meteorite and outgassing and lightning
    sp_names, sp_fluxes = dict_to_input(bombardment(bomb_rate[i], m_mass, prod_sp))
    # need: species, flux, vdep and name of txt file
    BC_bot_file = main_folder + 'BC_files/' + 'BC_bot_' + sim + '.txt'
    subprocess.check_call(['python', 'gen_BC.py', sp_names, sp_fluxes, vdep, BC_bot_file])
    # then use the new BC file in cfg file with output name
    BC_bot_change = 'bot_BC_flux_file,'+BC_bot_file+',str'
    out_file = sim + '.vul'
    out_change = 'out_name,'+out_file+',str'
    new_cfg = sim_folder + '/vulcan_cfg.py'
    subprocess.check_call(['python', 'gen_cfg.py', new_cfg, BC_bot_change, out_change])
    # then change to simulation folder and put symlinks in there to avoid copyying and make importing possible
    #subprocess.check_call(['cp', '-p', 'build_atm.py', 'chem_funs.py', 'op.py', 'phy_const.py', 'store.py', 'vulcan.py', sim_folder])
    wd = os.getcwd()
    os.chdir(sim_folder)
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/build_atm.py', 'build_atm.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/make_chem_funs.py', 'make_chem_funs.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/op.py', 'op.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/phy_const.py', 'phy_const.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/store.py', 'store.py'])
    subprocess.check_call(['cp', '-p', '/home/s2555875/VULCAN-2/vulcan.py', 'vulcan.py'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/thermo', 'thermo'])
    subprocess.check_call(['ln', '-s', '/home/s2555875/VULCAN-2/atm', 'atm'])
    
    # then run vulcan.py
    subprocess.check_call(['python', 'vulcan.py'])
    # then exit simulation folder and delete it
    os.chdir(wd)
    subprocess.check_call(['rm', '-rf', sim_folder])

# TO DO:
#   need something to prepare atmosphere and find a good place to store that, chemistry will be the same
# DONE  similar for boundary conditions
# DONE  make a generalised, probably with inputs, script for initial mixing ratios
# (AlMOST) DONE  write general vulcan_cfg generator -> NEED generic file to start from that has save structure
# DONE  do I need to copy vulcan.py for each function or can I call the same file? -> can call the same as they run as different processes with their own memory, variables, etc.
#   after all these can the parallelisation run
# DONE  plus figure out the output structure: how to name and what groupings should I save them in
# DONE  so far every bit is in a different script, maybe can be as functions here, THINK about this -> less command line python stuff if they're here as functions...