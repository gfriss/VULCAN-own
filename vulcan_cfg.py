# ============================================================================= 
# Configuration file of VULCAN:  
# ============================================================================= 

# ====== Setting up the elements included in the network ======
atom_list = ['H', 'O', 'C', 'N']
use_lowT_limit_rates = False

# ====== Setting up paths and filenames for the input and output files  ======
# input:
#network = 'thermo/CRAHCNO_V3.txt'
network = 'thermo/NCHO_full_photo_network.txt'
use_lowT_limit_rates = False
gibbs_text = '/home/s2555875/VULCAN-2/thermo/gibbs_text.txt' # (all the nasa9 files must be placed in the folder: thermo/NASA9/)
cross_folder = '/home/s2555875/VULCAN-2/thermo/photo_cross/'
#com_file = 'thermo/all_compose_crahcno.txt'
com_file = 'thermo/all_compose.txt'
atm_file = '/scratch/s2555875/TP_files/archean.txt' # TP and Kzz (optional) file
sflux_file = '/scratch/s2555875/stellar_flux/early_sun.txt' # This is the flux density at the stellar surface
top_BC_flux_file = '/home/s2555875/VULCAN-2/atm/BC_top_Earth.txt' # the file for the top boundary conditions
bot_BC_flux_file = '/home/s2555875/VULCAN-2/atm/BC_bot_Pearce_B.txt' # the file for the lower boundary conditions
vul_ini = '/home/s2555875/VULCAN-2/atm/mixing_table_archean.txt' # the file to initialize the abundances for ini_mix = 'vulcan_ini'
#vul_ini = '/scratch/s2555875/output/archean.vul'
# output:
output_dir = '/home/s2555875/scratch/output/'
plot_dir = '/home/s2555875/scratch/plot/'
movie_dir = '/home/s2555875/scratch/plot/movie/'
out_name = 'archean_ncho_nowash_test.vul' # output file name


# ====== Setting up the elemental abundance ======
use_solar = False # True: using the solar abundance from Table 10. K.Lodders 2009; False: using the customized elemental abundance. 
# customized elemental abundance (only read when use_solar = False)
O_H = 1.9185E+3#6.0618E-4 *(0.85) #*(0.793)  
C_H = 9.5925E+2#2.7761E-4  
N_H = 2.1316E+2#8.1853E-5
S_H = 2.1527E-7#1.3183E-5
He_H = 0.09692
ini_mix = 'table' # Options: 'EQ', 'const_mix', 'vulcan_ini', 'table' (for 'vulcan_ini, the T-P grids have to be exactly the same)

# Initialsing uniform (constant with pressure) mixing ratios (only reads when ini_mix = const_mix)
#const_mix = {'N2':0.78, 'O2':0.20, 'H2O':1e-6,  'CO2':4E-4, 'Ar':9.34e-3, 'SO2': 2e-10} 
#const_mix = {'H2':0.9, 'N2':0.1, 'H2O':1e-4, 'CH4': 2e-6, 'Ar':9.34e-3, 'SO2': 2e-10} # Ar and SO2 left from original
const_mix = {'N2':8.929e-01, 'CO2':0.1, 'H2O':4.112e-03, 'CH4':3e-3, 'O2':1e-7}

# ====== Setting up photochemistry ======
use_photo = True
# astronomy input
r_star = 1 # stellar radius in solar radius
Rp = 6.3781e8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 1 # planet-star distance in A.U.
sl_angle = 58 /180.*3.14159 # the zenith angle of the star in degree (usually 58 deg for the dayside average)
f_diurnal = 0.5 # to account for the diurnal average of solar flux (i.e. 0.5 for Earth; 1 for tidally-locked planets) 
scat_sp = ['N2', 'CO2'] # the bulk gases that contribute to Rayleigh scattering
T_cross_sp = ['CO2','H2O','NH3'] # warning: slower start! available atm: 'CO2','H2O','NH3', 'SH','H2S','SO2', 'S2', 'COS', 'CS2'

edd = 0.5 # the Eddington coefficient 
dbin1 = 0.1  # the uniform bin width < dbin_12trans (nm)
dbin2 = 2.   # the uniform bin width > dbin_12trans (nm)
dbin_12trans = 240. # the wavelength switching from dbin1 to dbin2 (nm)

# the frequency to update the actinic flux and optical depth
ini_update_photo_frq = 100
final_update_photo_frq = 5

# ====== Setting up ionchemistry ======
use_ion = False
if use_photo == False and use_ion == True:
    print ('Warning: use_ion = True but use_photo = False')
# photoionization needs to run together with photochemistry


# ====== Setting up parameters for the atmosphere ======
atm_base = 'N2' #Options: 'H2', 'N2', 'O2', 'CO2 -- the bulk gas of the atmosphere: changes the molecular diffsion, thermal diffusion factor, and settling velocity
rocky = True # for the surface gravity
nz = 120   # number of vertical layers
P_b = 1e6  # pressure at the bottom (dyne/cm^2)
P_t = 5e-2 # pressure at the top (dyne/cm^2)
use_Kzz = True
use_moldiff = True
use_vz = False
atm_type = 'file'  # Options: 'isothermal', 'analytical', 'file', or 'vulcan_ini' 'table'
Kzz_prof = 'file' # Options: 'const','file' or 'Pfunc' (Kzz increased with P^-0.4)
K_max = 1e5        # for Kzz_prof = 'Pfunc'
K_p_lev = 0.1      # for Kzz_prof = 'Pfunc'
vz_prof = 'const'  # Options: 'const' or 'file'
gs = 980.         # surface gravity (cm/s^2)  (HD189:2140  HD209:936)
Tiso = 1000 # only read when atm_type = 'isothermal'
# setting the parameters for the analytical T-P from (126)in Heng et al. 2014. Only reads when atm_type = 'analytical' 
# T_int, T_irr, ka_L, ka_S, beta_S, beta_L
para_warm = [120., 1500., 0.1, 0.02, 1., 1.]
para_anaTP = para_warm
const_Kzz = 1.E10 # (cm^2/s) Only reads when use_Kzz = True and Kzz_prof = 'const'
const_vz = 0 # (cm/s) Only reads when use_vz = True and vz_prof = 'const'

# frequency for updating dz and dzi due to change of mu
update_frq = 100 

# ====== Setting up the boundary conditions ======
# Boundary Conditions:
use_topflux = False
use_botflux = True
use_fix_sp_bot = {} #{"H2O":0.00894, 'CO2':4E-4, "H2O_l_s":0} #   0.0143 for 40% humidity 0.0033 for 20% humidity in US standard 1967 # fixed mixing ratios at the lower boundary
diff_esc = ['H2', 'H'] # species for diffusion-limit escape at TOA
max_flux = 1e13  # upper limit for the diffusion-limit fluxes

# ====== Reactions to be switched off  ======
remove_list = [] # in pairs e.g. [1,2]

# == Condensation ======
use_condense = True
use_settling = True
use_rainout = True
use_relax = ['H2O']
humidity = 0.25 # only for water
r_p = {'H2O_l_s': 0.01}  # particle radius in cm (1e-4 = 1 micron)
rho_p = {'H2O_l_s': 0.9}  # particle density in g cm^-3
start_conden_time = 0
condense_sp = ["H2O"]      
non_gas_sp = [ 'H2O_l_s']
rain_sp = ['H2O', 'HCN']
non_gas_rain_sp = ['H2O_rain', 'HCN_rain']
Henrys = 7.577 # M/atm
#fix_species = ['H2O','H2O_l_s'] # fixed the condensable species after condensation-evapoation EQ has reached  
fix_species = [] # no fixing due to rainout
fix_species_time = 5e8 # after this time to fix the condensable species
use_ini_cold_trap = False
stop_conden_time = 5e20 # making sure it does not happen
fix_species_from_coldtrap_lev = False

# ====== steady state check ======
st_factor = 0.5
conv_step = 500 

# ====== Setting up numerical parameters for the ODE solver ====== 
ode_solver = 'Ros2' # case sensitive
use_print_prog = True
use_print_delta = False
print_prog_num = 500  # print the progress every x steps 
dttry = 1.E-10
trun_min = 1e2
runtime = 1.E22
dt_min = 1.E-14
dt_max = runtime*1e-5
dt_var_max = 2.
dt_var_min = 0.5
count_min = int(2E2) #120 # already built in minimum countto check
count_max = int(1E5) #int(3E4)
atol = 1.E-5 # Try decreasing this if the solutions are not stable
mtol = 1.E-22
mtol_conv = 1.E-16
pos_cut = 0
nega_cut = -1.
loss_eps = 1e12 # for using BC
yconv_cri = 0.01 # for checking steady-state
slope_cri = 1.e-4
yconv_min = 0.1
flux_cri = 0.1
flux_atol = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)
ignore_layers = False
nl_ignore = 5
longdy_ignore = 1.

# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.05              # relative tolerence for adjusting the stepsize 
post_conden_rtol = 0.01 # switched to this value after fix_species_time

# ====== Setting up for ouwtput and plotting ======
# plotting:
plot_TP = False
use_live_plot = False
use_live_flux = False
use_plot_end = False
use_plot_evo = False
use_save_movie = False
use_flux_movie = False
plot_height = False
use_PIL = True 
live_plot_frq = 10
save_movie_rate = live_plot_frq
y_time_freq = 1  #  storing data for every 'y_time_freq' step
plot_spec = ['H2O', 'CO2', 'CO', 'CH4', 'NH3' ,'HCN', 'C2H2']  
# output:
output_humanread = False
use_shark = False
save_evolution = True   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
save_evo_frq = 100
save_if_converged = True