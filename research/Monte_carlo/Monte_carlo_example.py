#%% ----------------------------------------------------------------------------
# A. Hockin, January 2023
# KWR 403230-003
# Pipe permeation calculator
# With Martin vd Schans, Bram Hillebrand, Lennart Brokx
# ------------------------------------------------------------------------------

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from datetime import timedelta
from scipy.optimize import minimize
import itertools
from project_path import file_path
import matplotlib.pyplot as plt
import pickle
import os
from tqdm import tqdm #tqdm gives a progress bar for the simultation
from statistics import NormalDist
import random

from pipepermcalc.pipe import * 
from pipepermcalc.segment import * 

def check_create_folders(folder_name):
    # Check if folder exists, if not create it
    MYDIR = (folder_name)
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        # print("created folder : ", MYDIR)
    # else:
    #     print(MYDIR, "folder already exists.")
    save_results_to = MYDIR
    return MYDIR

def save_df_pickle(filename, df, foldername=None):
    #Save to pickle files 
    filename = filename+'.pickle'
    if foldername:
        outfile = open(foldername+'/'+filename,'wb')
    else: outfile = open(filename,'wb')
    pickle.dump(df,outfile)
    outfile.close()

def load_pickle(filename, foldername=None):
    filename = filename+'.pickle'
    if foldername:
        infile = open(foldername+'/'+filename,'rb')
    else: infile = open(filename,'rb')

    df = pickle.load(infile)
    infile.close
    return df

#%%

seg1 = Segment(name='seg1',
            material='PE40',
            length=25, #here set the length
            inner_diameter=0.0196,
            wall_thickness=0.0027)

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_conditions(concentration_soil = 15, #here set the gw conc
                chemical_name="Benzeen", 
                temperature_groundwater=12,
                flow_rate=1.5, 
                suppress_print=True)
seg1.partitioning_enthalpie
#%%
# Overview of steps for Monte-Carlo simulations

#Parameters to vary:
    # size of plume (length of segment) (excel file) -> Amitosh/Auliato discuss more complicated version later
    # concentration of plume (excel file)
    # partitioning coefficient, K_ref, +/- 0.5 st. dev, @Martin to give the real values later
    # diffusion coefficient, K_ref, +/- 0.5 st. dev, @Martin to give the real values later
    # Assessment factor = 1 or 3
    # Flow rate -> @Martin, contact dw companies for distribution 
    # Pipe thickness ? -> Amitosh/Aulia

# Steps
# import the data on plume concentration, contact length
# set range for diffusion coefficient, D_ref +/- XX ->  See table 5-1 KWR 2016.056
# set range for partitioning coefficient, K_ref +/- XX ->  See table 5-1 KWR 2016.056
# set range flow rate
# set range pipe thickness

# Normal distrbution from mean (mu) and standard deviation (sigma), p is the 
# random number from the random number generator between 0 and 1
# NormalDist(mu, sigma).inv_cdf(p)

#%%
# import the data on plume concentration
# df = pd.read_excel(module_path / 'research' / 'Monte_carlo' / '20190702 kans normoverschrijding.xlsx', 
#                     sheet_name='RGW_AH', header=[175], usecols = "A:E", ) 

# save_df_pickle(filename='monte-carlo_plume_concs', df= df, foldername='monte-carlo_output')
df = load_pickle(filename='monte-carlo_plume_concs', foldername='monte-carlo_output')
ext_values = list(df.ext_value)

# range lenth pipe 
# length_range =pd.read_excel(module_path / 'research' / 'Monte_carlo' / '20190702 kans normoverschrijding.xlsx', 
#                     sheet_name='RGW_AH', header=[0], usecols = "U", nrows=12) 
# save_df_pickle(filename='monte-carlo_lenths', df= length_range, foldername='monte-carlo_output')
length_range = load_pickle(filename='monte-carlo_lenths', foldername='monte-carlo_output')
length_values = list(length_range.contactlengte)

#%%
save_results_to = check_create_folders(folder_name='monte-carlo_output')

# Loop through the combinations, save the dw concentrations
dw_concs = []
soil_concs = []
log_Kpws = []
log_Dps = []
lengths = []
flow_rates = []
wall_thicknesses =[]

log_Dp_refs =[]
log_Kpw_refs =[]
f_Dconcs =[]
f_Kconcs =[]
f_Dtemps =[]
f_Ktemps =[]

# initialize the index parameters
tenth_perc_n_min_1 = 0
ninety_perc_n_min_1 = 0

# Random number seeded to always produce the same sequence of random numbers
random.seed(5) 

# Set number of simulations per round, set tolerance for checking simulation rounds
sims = range(1000)
tolerance = 0.01

# Choice of assessment factor, 3 or 1
assessment_factor = 3

while True:

    for lp in tqdm(sims):
        # Input variables
        # ---------------
        soil_conc = random.choice(ext_values)

        # ah_todo
        # Update with information from Amitosh/Aulia
        length = random.choice(length_values)

        #ah_todo
        # Update with information from Mirjam/DWC
        flow_rate = NormalDist(mu=0.5, sigma=0.1).inv_cdf(p=random.random())
        #ah_todo change this to be the 1,2,4 person households at 120 L per person per day

        # ah_todo
        # Update with information from Amitosh/Aulia
        wall_thickness = NormalDist(mu=0.0027, sigma=0.0001).inv_cdf(p=random.random())

        seg1 = Segment(name='seg1',
                        material='PE40',
                        length=length, #here set the length
                        inner_diameter=0.0196,
                        wall_thickness=wall_thickness)

        pipe1 = Pipe(segment_list=[seg1])

        pipe1.set_conditions(concentration_soil = soil_conc, #here set the gw conc
                            chemical_name="Benzeen", 
                            temperature_groundwater=12,
                            flow_rate=flow_rate, 
                            suppress_print=True)

        # Update the partitioning and diffusion coefficients
        # --------------------------------------------------
        #Sr = standard error of regression
        # Values from 20160703 Database labmetingen

        #Reference D, K values
        Sr_D = 0.19572320 #Table 5-5, KWR 2016.056, excel:'PermDbase' DM-25
        log_Dp_ref = NormalDist(mu=seg1.log_Dp_ref, #-11.54717333172 #
                                sigma=Sr_D).inv_cdf(p=random.random())
        
        Sr_K = 0.31397266 #Table 5-5, KWR 2016.056, excel:'PermDbase' AL-25
        log_Kpw_ref = NormalDist(mu=seg1.log_Kpw_ref, #1.6476099999999998 #
                                sigma=Sr_K).inv_cdf(p=random.random())
        
        # concentration corrections
        Sr_conc_D = 0.07662645 #excel:'CONC' AE-4
        f_Dconc = NormalDist(mu=seg1.f_Dconc, 
                            sigma=Sr_conc_D).inv_cdf(p=random.random())
        
        Sr_conc_K = 0.10106212 #excel:'CONC' W-4
        f_Kconc = NormalDist(mu=seg1.f_Kconc, 
                            sigma=Sr_conc_K).inv_cdf(p=random.random())

        # temperature corrections
        #@MartinvdS corrections on the act. engergy/enthalpie itself not the factor
        Sr_temp_D = 11.7958431 #Table 5-6, KWR 2016.056, excel:'TEMP' CO-125
        activattion_energy = NormalDist(mu=seg1.activattion_energy, 
                            sigma=Sr_temp_D).inv_cdf(p=random.random())
        
        f_Dtemp = (activattion_energy / (0.008314 * np.log(10)) 
                  * (1 / (25 + 273) - 1 / (pipe1.temperature_groundwater + 273)))

        Sr_temp_K = 13.2239059 #Table 5-6, KWR 2016.056, excel:'TEMP' CJ-125
        partitioning_enthalpie = NormalDist(mu=seg1.partitioning_enthalpie, 
                            sigma=Sr_temp_K).inv_cdf(p=random.random())
        
        f_Ktemp = (partitioning_enthalpie / (0.008314 * np.log(10)) 
                  * (1 / (25 + 273) - 1 / (pipe1.temperature_groundwater + 273)))

        # age corrections @MartinvdS include the age corrections?
        Sr_age_D = 0.17 #Eqn 21, KWR 2016.056, excel:?
        f_Dage = NormalDist(mu=0, 
                            sigma=Sr_age_D).inv_cdf(p=random.random())
        
        Sr_age_K = 0.05 #Eqn 22, KWR 2016.056, excel:?
        f_Kage = NormalDist(mu=0, 
                            sigma=Sr_age_K).inv_cdf(p=random.random())

        # Set the Kpw and Dp
        seg1.log_Kpw = log_Kpw_ref + f_Kconc + f_Ktemp + f_Kage
        seg1.log_Dp = log_Dp_ref + f_Dconc + f_Dtemp + f_Dage

        #set assessment_factor
        pipe1.ASSESSMENT_FACTOR_GROUNDWATER = assessment_factor

        pipe1.validate_input_parameters()
    
        # Calculate concentrations, can we do in one loop and store seperate peak/mean conc
        # dw_conc = pipe1.calculate_mean_dw_concentration()
        dw_conc = pipe1.calculate_peak_dw_concentration()

        dw_concs.append(dw_conc)
        soil_concs.append(pipe1.concentration_groundwater)  
        lengths.append(seg1.length)
        flow_rates.append(pipe1.flow_rate)
        wall_thicknesses.append(seg1.wall_thickness)
        log_Dp_refs.append(log_Dp_ref)
        log_Kpw_refs.append(log_Kpw_ref)
        f_Dconcs.append(f_Dconc)
        f_Kconcs.append(f_Kconc)
        f_Dtemps.append(f_Dtemp)
        f_Ktemps.append(f_Ktemp)
        log_Kpws.append(seg1.log_Kpw)  
        log_Dps.append(seg1.log_Dp)

    # check if the 10th and 90th percentile within tolerance, then stop the loop
    tenth_perc = np.percentile(dw_concs, 10)
    ninety_perc = np.percentile(dw_concs, 90)

    criteria_ten = abs(1 - tenth_perc / tenth_perc_n_min_1)
    criteria_nine = abs(1 - ninety_perc / ninety_perc_n_min_1)

    if (criteria_ten <= tolerance) and (criteria_nine <= tolerance):
        break
    elif len(dw_concs) > 100000: # break if the code takes too many simulations
        break
    else: 
        ninety_perc_n_min_1 = ninety_perc
        tenth_perc_n_min_1 = tenth_perc

    print('ninety_perc:', ninety_perc, 'tenth_per:c', tenth_perc)

# put the data into a df, sort by the dw_conc and save
data = zip(dw_concs, soil_concs, log_Kpws, log_Dps, lengths, flow_rates, wall_thicknesses, 
           log_Dp_refs, log_Kpw_refs, f_Dconcs, f_Kconcs, f_Dtemps, f_Ktemps)
df = pd.DataFrame(data,  columns = ['dw_concs', 'soil_conc', 'Kpw', 'Dp', 'Length', 'flow_rates', 'wall_thicknesses',
                                    'log_Dp_refs', 'log_Kpw_refs','f_Dconcs', 'f_Kconcs', 'f_Dtemps', 'f_Ktemps'])
df.sort_values(by=['dw_concs'], inplace=True)
df.reset_index(inplace=True)

save_df_pickle(filename='example_monte-carlo_loop_df', df= df, foldername='monte-carlo_output')

dw_norm = 0.001 # Benzene drinking water norm, g/m3
print('Exceedences:', round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 3), '%, total sims:', len(df) )

#%%
save_results_to = check_create_folders(folder_name='figures')

df = load_pickle(filename='example_monte-carlo_loop_df', foldername='monte-carlo_output')
fig = plt.figure(figsize=[10, 5])

plt.plot(df.dw_concs, df.index/len(df), ) 
plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
plt.xlabel('Mean DW concentration (g/m3)')
plt.ylabel('Probability')
plt.title('Exceedences:'+str(round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 3))+ '%, total sims:'+ str(len(df)) )
plt.xscale('log')
plt.legend()
plt.savefig(save_results_to +'/example_peak_monte-carlo_simulation.png', dpi=300, bbox_inches='tight')

# %%
# Creating a series of data of in range of 1-10000.
lps = np.arange(start=0, stop=10000, step=1)
cdf_K = []
cdf_D = []

# calculate the cdf for the different parameters...
for lp in lps:

    Sr_K = 0.31397266
    log_Kpw_ref = NormalDist(mu=seg1.log_Kpw_ref, #1.6476099999999998 #
                             sigma=Sr_K).inv_cdf(p=random.random())

    Sr_D = 0.19572320 #Sr = standard error of regression
    log_Dp_ref = NormalDist(mu=seg1.log_Dp_ref, #-11.54717333172 #
                            sigma=Sr_D).inv_cdf(p=random.random())


    cdf_K.append(log_Kpw_ref )
    cdf_D.append(log_Dp_ref )

cdf_K.sort()
cdf_D.sort()

#Plotting the Results
# plt.plot(cdf_K, lps, color = 'red')
plt.plot(cdf_D, lps, color = 'blue')

plt.xlabel('log Kpw/Dp')
plt.ylabel('Cumulative Density')

np.mean(cdf_K), np.std(cdf_K), np.mean(cdf_D), np.std(cdf_D)
# %%
