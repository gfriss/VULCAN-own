#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

import pickle
#%%
def rainout(dat, rain_spec = 'HCN_rain', g_per_mol = 27):
    rain_rate = np.sum(dat['variable']['y_rain'][rain_spec][:-1] * dat['atm']['dzi']) / dat['variable']['dt'] # 1/cm2s
    rain_rate = rain_rate * 2.259e-13 # mol/m2yr
    return rain_rate * (g_per_mol/1000.) # kg/m2yr

nsim = 15 # hardcoded for now, change later...
out_folder = '/scratch/s2555875/output/'
bc_folder= '/scratch/s2555875/BC_files/'
plot_folder = '/scratch/s2555875/plot/'
hcn_influx = [] # list to store the results
prec = []
#%%
for i in range(nsim):
    sim = 'sim_' # setting up the names of files
    if i < 10:
        sim += '0' + str(i)
    else:
        sim += str(i)
    sim += '_onlyH2.vul'
    
    with open(out_folder+sim, 'rb') as handle: # reading in files
        data = pickle.load(handle)
        hcn_influx.append(rainout(data))
        prec.append(rainout(data, rain_spec = 'H2O_rain', g_per_mol = 18))


# post simulations...
nsim_extra = 8
hcn_influx_extra = []
prec_extra = []
for i in range(nsim_extra):
    sim = 'sim_' # setting up the names of files
    if i < 10:
        sim += '00' + str(i)
    else:
        sim += str(i)
    sim += '_onlyH2.vul'
    
    with open(out_folder+sim, 'rb') as handle: # reading in files
        data = pickle.load(handle)
        hcn_influx_extra.append(rainout(data))
        prec_extra.append(rainout(data, rain_spec = 'H2O_rain', g_per_mol = 18))


#hcn_influx = [hcn_influx[0]] + hcn_influx_extra[::4] + hcn_influx[1:][::6] # to have bombardment rates in order
hcn_influx_new = [hcn_influx[0]] + [hcn_influx_extra[3]] + [hcn_influx[2]]
hcn_influx = np.array(hcn_influx_new)
prec = [prec[0]] + prec_extra + prec[1:]
prec = np.array(prec)
#%%
sims = [out_folder + 'sim_07_onlyH2.vul', out_folder + 'sim_07_CtoO.vul']
hcn_influx.append(3.4e-12)
for s in sims:
    with open(s, 'rb') as handle:
        sim_data = pickle.load(handle)
    hcn_influx.append(rainout(sim_data))
hcn_influx = np.array(hcn_influx)
lab = ['Pearce et al. (2022)', 'Mass delivery rate = 5.15e+24 g/Gyr', 'C/O ratio = 0.664']
#%%
plt.clf()

#Variable Declarations
T_celsius = 327.7 - 273.15
T = 327.7

w_i = 60.7e-9
m_dot_I = 6e8
f_s = 0.32
r = 40.
rho = 2185.
r_p = 1.
A_p = math.pi*r_p**2
r_g = 500.
tau_d_1cm = 4.9e-3
tau_d_5cm = 0.12
tau_d_10cm = 0.48
R_plus = 6371000
gamma = 31557600
rho_w = 1000.
F = 0.4

tau_s = 1.
P = 3.5
delta = 0.5
sp = 0.3
min_water = 0.001 # 1mm
S = 0.95 #Seepage rate 0.36525
Phi = 1e-4
lambda_uv = 225e-9
h = 6.626e-34
c = 2.9979e8
N_A = 6.022e23 #Avogadro's number



#initialising species (going to be dict in dict)
species = ['HCN', 'adenine', 'guanine', 'uracil', 'cytosine', 'thymine', '2amino', 'ribose', 'formaldehyde', 'xanthine', 'hypoxanthine']

#and mu values in kg/mol
# 'HCN', 'adenine', 'guanine', 'uracil', 'cytosine', 'thymine', '2amino', 'ribose', 'formaldehyde', 'xanthine', 'hypoxanthine'
MU = [0.0270253, 0.13513, 0.15113, 0.1120868, 0.1111, 0.1261133, 0.084077, 0.15013, 0.030031, 0.15211, 0.1361115]

#and densities to calc d 
rho_species = [687, 1470, 2200, 1320, 1550, 1230, 800, 1200, 815, 1600, 2000]
def calc_d(mu_val, rho_val):
    ''' Calculates molecular density(?)'''
    return 2*(3*mu_val/(4*math.pi*N_A*rho_val))**(1./3)

# photo dissociation rates (kg/yr/m^2)
def photodis_rate(muval):
    ''' Calculates the photodissociation rate for a given molecule'''
    return ((Phi*F*lambda_uv*gamma*muval)/(h*c*N_A))

# k values, eq part for non-nucleobase is missing so k is set to 0 for such species
k = [0., 10**(-5902/T + 8.15), 10**(-6330/T + 9.40), 10**(-7649/T + 11.76), 10**(-5620/T + 8.69), 10**(-7709/T + 11.24), 0., 0., 0., 10**(-6230/T + 9.42), 10**(-5270/T + 7.95)]


def calc_m_dot(spec, mu_val, mu_val_hcn):
    ''' Calculates the mass increase due to the incoming mass flux of rain out of HCN and formaldehyde
        in the reducing model'''
    if spec == 'HCN':
        return hcn_influx * 4 * math.pi * R_plus**2
    #elif spec == 'formaldehyde':
    #    return H2CO_mass_influx_red * 4 *math.pi * R_plus**2
    else:
        return hcn_influx * 4 * math.pi * R_plus**2 * mu_val/mu_val_hcn
    

# initialising dict in dict to store the all the info
constants_and_rates = {}
for spec,mu_spec,rho_spec,k_spec in zip(species,MU,rho_species,k):
    constants_and_rates[spec] = {'mu': mu_spec,
                       'd': calc_d(mu_spec, rho_spec),
                       'rho': rho_spec,
                       'M_uv_dot': photodis_rate(mu_spec),
                       'm_dot': calc_m_dot(spec, mu_spec, MU[0]),
                       'k': k_spec}
    

#Fraction of surviving organics during entry
f_s = {'IDP': 0.06, 'Met': 0.32}

E = S-0.12 + 0.06*T_celsius
tmax = 8 #years
level = 16

nt = (2**level) + 1 #Choosing nt to have twice as many grid points as nx
    
# Array for plotting
t = np.linspace(0,tmax,nt)

# Calculate delta_t
delta_t = t[2] - t[1]

#Constant seepage mass per year
m_seepage_rate = math.pi*rho_w*r_p**2*S

m_i0 = (4./3)*w_i*f_s['Met']*r**3*rho*A_p/r_g**2

pause_Met = 0

#initialising, original mixing a lot oxidising for some reason and last two nucleotides...now included
values = {}
names = ['L_IDP', 'L_Met', 'm_IDP', 'm_Met', 'm_IDP_A', 'm_Met_A', 'C_IDP', 'C_Met'] +\
['m_' + spec for spec in species] + ['C_' + spec for spec in species]
for name in names:
    if name not in ['L_IDP', 'L_Met', 'm_IDP', 'm_Met', 'm_IDP_A', 'm_Met_A', 'C_IDP', 'C_Met']:
        values[name] = np.zeros(shape=(len(hcn_influx), nt))
    else:
        values[name] = np.zeros(shape=nt)
values['L_IDP'][0] = r_p - min_water
values['L_Met'][0] = r_p - min_water
values['m_IDP'][0] = math.pi*rho_w*r_p**2*(min_water) #og was r_p-(r_p-min_water)....
values['m_Met'][0] = math.pi*rho_w*r_p**2*(min_water)
    

def precipitation(iteration):
    ''' Calculates the precipitation rate that is sinusoidal to represent the seasonal cycle for the current iteration.
        P (mean precipitation rate) and sp (seasonal phase shift) are model dependent.'''
    return (delta_t*P)*(1 + delta*np.sin(2*math.pi*(t[iteration] - sp)/tau_s))

def water_decrease(impactor, iteration):
    ''' Calculates the rate of water decrese by adding up evaporation and seepage then subtraction the precipitation.
        This is also done for current iteration steps and all precipitation models.'''
    return E*delta_t + values['L_'+impactor][iteration] - precipitation(iteration)

def nucleobase_outflow(iteration, pause, kval, model):
    ''' Calculates the nucleobase outflow rate from meteorites for current iteration step. This only happens when the pond is wet, hence the introduction of the variable pause.
        In case of IDP this is set to zero as outflow from small IDPs are basically instantenaous compared to the simulation time.'''
    if model == 'IDP':
        return 0.
    elif model == 'Met':
        return delta_t * m_i0 * np.e**(-t[iteration-pause]*(gamma*kval + (1./tau_d_1cm))) / tau_d_1cm
    
def water_mass(impactor, iteration):
    ''' Calculates the current water mass.'''
    return math.pi*rho_w*r_p**2*(r_p-values['L_'+impactor][iteration])

def mass_increase(mdot, model = 'IDP', wi = 1., fs = 1.):
    ''' Calculates the mass increase using the timestep (delta_t) and unit pond surface area (A_p/(4pi*R_plus**2))
        This term is not present when working with Meteorites. For generalisation porpuses, though, it appaers there but is set to 0.'''
    if model == 'Met':
        return 0
    elif model == 'IDP':
        return (delta_t * wi * mdot * fs * A_p) / (4*np.pi * R_plus**2)
    
def photo_destruction(muvdot, mass, RHO, den):
    ''' Calculates the photo destruction rate depending on the amount of nucleobes:
        if there's enough nucleobase to cover the base of the pond then we multiply by A_p;
        if there isn't, we multiply by the cross-sectional area'''
    if mass / (RHO*den) < A_p:
        return -delta_t * muvdot * mass / (RHO*den)
    else:
        return -delta_t * muvdot * A_p

def decomposition(K_val, mass):
    ''' Calculates the decomposition/hydrolisis rate using the rate constant (k) and the constant gamma'''
    return -delta_t * gamma * K_val * mass

def seepage(mass, mass_Met_or_IDP):
    ''' Calculates the seepage rate for each specias using a constant overall seepage rate'''
    return -delta_t*mass*m_seepage_rate/mass_Met_or_IDP


# Solve ODE numerically
# Biomolecule evolution from meteorites, IDPs, and aqueous production from atmospheric precursors
for n in range(0,nt-1):
    # first cycle through the Meteorite and IDP equations
    for impact in ['Met', 'IDP']:
        values['L_'+impact][n+1] = water_decrease(impact, n)
        if (values['L_'+impact][n+1] < 0):
            values['L_'+impact][n+1] = 0
        if (values['L_'+impact][n+1] >= (r_p - min_water)):
            values['L_'+impact][n+1] = r_p - min_water
            values['m_'+impact+'_A'][n+1] = values['m_'+impact+'_A'][n] + mass_increase(m_dot_I, model = impact, wi = w_i, fs = f_s[impact]) + photo_destruction(constants_and_rates['adenine']['M_uv_dot'], values['m_'+impact+'_A'][n], constants_and_rates['adenine']['rho'], constants_and_rates['adenine']['d'])
            if impact == 'Met':
                pause_Met += 1
        else:
            values['m_'+impact+'_A'][n+1] = values['m_'+impact+'_A'][n] + mass_increase(m_dot_I, model = impact, wi = w_i, fs = f_s[impact]) + nucleobase_outflow(n, pause_Met, constants_and_rates['adenine']['k'], model = impact) + decomposition(constants_and_rates['adenine']['k'], values['m_'+impact+'_A'][n]) + seepage(values['m_'+impact+'_A'][n], values['m_'+impact][n])
        if values['m_'+impact+'_A'][n+1] < 0:
            values['m_'+impact+'_A'][n+1] = 0
        values['m_'+impact][n+1] = water_mass(impact, n+1)
        # update concentration
        values['C_'+impact][n+1] = values['m_'+impact+'_A'][n+1] / values['m_'+impact][n+1]    
    # then the aqueous equations
    for name in names:
        for i in range(len(hcn_influx)):
            # we only work with masses in the pond, not from Met or IDP here (= aqueous production)
            if name[0] == 'm' and any([x in name for x in ['Met', 'IDP']]) == False:
                #print(name)
                spec_name = name.split('_')[1] # the format is m_spec
                spec_mdot = constants_and_rates[spec_name]['m_dot'][i] # this is a list with values from all influxes
                
                if (values['L_IDP'][n+1] >= (r_p - min_water)):
                    values['L_IDP'][n+1] = r_p - min_water
                    values[name][i, n+1] = values[name][i, n] + mass_increase(spec_mdot) + photo_destruction(constants_and_rates[spec_name]['M_uv_dot'], values[name][i, n], constants_and_rates[spec_name]['rho'], constants_and_rates[spec_name]['d'])
                else:
                    values[name][i, n+1] = values[name][i, n] + mass_increase(spec_mdot) + decomposition(constants_and_rates[spec_name]['k'], values[name][i, n]) + seepage(values[name][i, n], values['m_IDP'][n])
                if values[name][i, n+1] < 0:
                    values[name][i, n+1] = 0
                # update concentrations (biomolecule mass / water mass)
                conc_name = 'C' + name[1:]
                values[conc_name][i, n+1] = values[name][i, n+1] / values['m_IDP'][n+1]
            


#Conversion from molar to mass mixing ratios
def molar2mass(x):
    return x * 1.e3 * constants_and_rates['adenine']['mu']

def mass2molar(x):
    return x / 1.e3 / constants_and_rates['adenine']['mu']

#Experimental yields
Adenine_lower = 0.005
Adenine_upper = 0.18
Guanine_lower = 6.7e-5
Guanine_upper = 0.2
Cytosine = 0.036
Uracil_lower = 1.7e-5
Uracil_upper = 0.018
Thymine = 0.012
Two_Amino_oxazole = 0.0011
Ribose = 3.6e-4
Formaldehyde = 0.036
#%%
f, ax1 = plt.subplots(figsize = (12,10))
#f, (ax1, ax2) = plt.subplots(1, 2, figsize=(28,10))
        
for i in range(len(hcn_influx)):
    p1 = ax1.fill_between(t, values['C_adenine'][i]*Adenine_lower*1e6/constants_and_rates['adenine']['mu'], values['C_adenine'][i]*Adenine_upper*1e6/constants_and_rates['adenine']['mu'], linestyle='-', lw=3.5, alpha=.40, label = lab[i])#r'$k_{HCN rain} = $' + '{:.2e}'.format(hcn_influx[i]) + r' kg m$^{-2}$ yr$^{-1}$')
              
secax = ax1.secondary_yaxis('right', functions=(molar2mass, mass2molar))

secax.set_ylabel("Adenine Mass Fraction (ppb)",fontsize=22)

#ax1.set_ylim(5e-8,1e-1)

ax1.set_yscale('log')
ax1.legend(ncol = 2, fontsize = '16', loc = 'lower center')

  
ax1.set_xlabel('Time (yr)', fontsize=22)
ax1.set_ylabel('Adenine Molar Concentration ($\mu$M)', fontsize=22)
for tick in ax1.xaxis.get_major_ticks():
    tick.label1.set_fontsize(22) 
for tick in ax1.yaxis.get_major_ticks():
    tick.label1.set_fontsize(22)

ax1.tick_params(which='both', direction='out', length=6, width=2)
secax.tick_params(which='both', direction='out', length=6, width=2)
secax.tick_params(axis='both',labelsize=22)

f.savefig('/scratch/s2555875/plot/adenine_pres.pdf', bbox_inches = 'tight')

# %%
