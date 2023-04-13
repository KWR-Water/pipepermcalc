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
# Overview of steps for Monte-Carlo simulations

#Parameters to vary:
    # size of plume (length of segment) (excel file)
    # concentration of plume (excel file)
    # partitioning coefficient, K_ref, +/- st. dev, @Martin to give the real values later
    # diffusion coefficient, K_ref, +/- st. dev, @Martin to give the real values later
    # Assessment factor = 1 or 3
    # Flow rate -> @Martin, contact dw companies for distribution 
    # Pipe thickness ? -> Amitosh/Aulia

# Normal distrbution from mean (mu) and standard deviation (sigma), p is the 
# random number from the random number generator between 0 and 1
# NormalDist(mu, sigma).inv_cdf(p)

#%%
# import the data on plume concentration
df = pd.read_excel('20190702 kans normoverschrijding - Copy.xlsx', 
                    sheet_name='`Brondata anoniem', header=[3], usecols = "P", ) 
save_df_pickle(filename='monte-carlo_plume_concs', df= df, foldername='monte-carlo_output')

#%%
# range lenth plume
# plume_length_range =pd.read_excel('20190702 kans normoverschrijding.xlsx', 
#                     sheet_name='RGW_AH', header=[0], usecols = "U", nrows=12) 
# save_df_pickle(filename='monte-carlo_lenths', df= plume_length_range, foldername='monte-carlo_output')
plume_length_range = load_pickle(filename='monte-carlo_lenths', foldername='monte-carlo_output')
plume_length_values = list(plume_length_range.contactlengte)

#%%
# import data on pipe dimentsion and concentrations
df = load_pickle(filename='monte-carlo_plume_concs', foldername='monte-carlo_output')
plume_concs = list(df['Extrapolated values'])

# Import the pipe information - inner diameter, thickness, length
df_PE40 = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40')

# range lenth pipe
df_PE40_Lengte_GIS = df_PE40.sort_values(by='Lengte_GIS')
length_values = list(df_PE40_Lengte_GIS.Lengte_GIS)

# range inner diameter
df_PE40_inner_diam = df_PE40[df_PE40['Binnendiam'] > 0]  
df_PE40_inner_diam = df_PE40_inner_diam.sort_values(by='Binnendiam')
inner_diam_values = list(df_PE40_inner_diam.Binnendiam)

# dictionary relating specific inner diameters to specific wall thicknesses    
wall_thickness_dict = {12.39999962: 1.8,
                        15.60000038: 2.2, 
                        19.60000038: 2.7,
                        25: 3.5,
                        31.39999962: 4.3,
                        39.20000076: 5.4,
                        49.40000153: 6.8,}

#%%
save_results_to = check_create_folders(folder_name='monte-carlo_output')

def run_Monte_Carlo_simulation (plume_concs, 
                                length_values,
                                inner_diam_values,
                                wall_thickness_dict,
                                assessment_factor,
                                calculate_peak=False, 
                                calculate_mean=False,):

    # Loop through the combinations, save the dw concentrations
    dw_concs = []
    soil_concs = []
    gw_concs = []
    log_Kpws = []
    log_Dps = []
    lengths = []
    flow_rates = []
    inner_diameters = []
    wall_thicknesses =[]

    log_Dp_refs =[]
    log_Kpw_refs =[]
    f_Dconcs =[]
    f_Kconcs =[]
    activattion_energies = []
    partitioning_enthalpies = []
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

    while True:

        for lp in tqdm(sims):
            # Input variables
            # ---------------
            # Soil Concentration, NBNL data
            soil_conc = random.choice(plume_concs)

            # Length, PWN data given in meters
            length = random.choice(length_values)

            # Inner diameter, PWN data, given in mm (convert to m)
            inner_diameter_mm = random.choice(inner_diam_values) 
            inner_diameter = inner_diameter_mm / 1000

            # matching wall thicnkess from inner diameter, PWN data, given in mm (convert to m)
            wall_thickness_mm = (wall_thickness_dict[inner_diameter_mm])
            wall_thickness = wall_thickness_mm / 1000

            #ah_todo
            # Update with information from Mirjam/DWC ? @MartinvdS
            flow_rate = NormalDist(mu=0.5, sigma=0.1).inv_cdf(p=random.random())
            #ah_todo change this to be the 1,2,4 person households at 120 L per person per day?

            # Create pipe and set conditions
            # ------------------------------
            seg1 = Segment(name='seg1',
                            material='PE40',
                            length=length,
                            inner_diameter=inner_diameter,
                            wall_thickness=wall_thickness)

            pipe1 = Pipe(segment_list=[seg1])

            pipe1.set_conditions(concentration_soil = soil_conc, #here set the gw conc
                                chemical_name="Benzeen", 
                                temperature_groundwater=12,
                                flow_rate=flow_rate, 
                                suppress_print=True)

            # Update the partitioning and diffusion coefficients
            # --------------------------------------------------
            # Sr = standard error of regression
            # Values from 20160703 Database labmetingen excel file

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
            Sr_age_D = 0.17 #Eqn 21, KWR 2016.056
            f_Dage = NormalDist(mu=0, 
                                sigma=Sr_age_D).inv_cdf(p=random.random())
            
            Sr_age_K = 0.05 #Eqn 22, KWR 2016.056
            f_Kage = NormalDist(mu=0, 
                                sigma=Sr_age_K).inv_cdf(p=random.random())

            # Set the Kpw and Dp
            seg1.log_Kpw = log_Kpw_ref + f_Kconc + f_Ktemp + f_Kage
            seg1.log_Dp = log_Dp_ref + f_Dconc + f_Dtemp + f_Dage

            #set assessment_factor
            pipe1.ASSESSMENT_FACTOR_GROUNDWATER = assessment_factor

            pipe1.validate_input_parameters()
        
            # Calculate concentrations, can we do in one loop and store seperate peak/mean conc
            if calculate_peak:
                dw_conc = pipe1.calculate_peak_dw_concentration()
            elif calculate_mean:
                dw_conc = pipe1.calculate_mean_dw_concentration()
                    
            # Save the calculation information
            dw_concs.append(dw_conc)
            soil_concs.append(pipe1.concentration_soil)  
            gw_concs.append(pipe1.concentration_groundwater)  
            lengths.append(seg1.length)
            flow_rates.append(pipe1.flow_rate)
            inner_diameters.append(inner_diameter)
            wall_thicknesses.append(seg1.wall_thickness)
            log_Dp_refs.append(log_Dp_ref)
            log_Kpw_refs.append(log_Kpw_ref)
            f_Dconcs.append(f_Dconc)
            f_Kconcs.append(f_Kconc)
            f_Dtemps.append(f_Dtemp)
            f_Ktemps.append(f_Ktemp)
            activattion_energies.append(activattion_energy)
            partitioning_enthalpies.append(partitioning_enthalpie)
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
    data = zip(dw_concs, soil_concs, gw_concs, log_Kpws, log_Dps, lengths, 
               inner_diameters, flow_rates, wall_thicknesses, 
            log_Dp_refs, log_Kpw_refs, f_Dconcs, f_Kconcs, f_Dtemps, f_Ktemps, 
            activattion_energies, partitioning_enthalpies)
    df = pd.DataFrame(data,  columns = ['dw_concs', 'soil_conc', 'gw_concs', 'Kpw', 'Dp', 'Length', 'inner_diameter', 
                                        'flow_rates', 'wall_thicknesses',
                                        'log_Dp_refs', 'log_Kpw_refs','f_Dconcs', 'f_Kconcs', 'f_Dtemps', 'f_Ktemps',
                                          'activattion_energies', 'partitioning_enthalpies'])
    df.sort_values(by=['dw_concs'], inplace=True)
    df.reset_index(inplace=True)

    return df
#%%
df_peak = run_Monte_Carlo_simulation (plume_concs=plume_concs, 
                                length_values=length_values,
                                inner_diam_values=inner_diam_values,
                                wall_thickness_dict=wall_thickness_dict,
                                assessment_factor=3,
                                calculate_peak=True,)

save_df_pickle(filename='example_peak_monte-carlo_loop_df', df= df_peak, foldername='monte-carlo_output')

dw_norm = 0.001 # Benzene drinking water norm, g/m3
print('Peak exceedences:', round(len(df_peak.loc[df_peak.dw_concs > dw_norm]) / len(df_peak)*100, 3), '%, total sims:', len(df_peak) )


df_mean = run_Monte_Carlo_simulation (plume_concs=plume_concs, 
                                length_values=length_values,
                                inner_diam_values=inner_diam_values,
                                wall_thickness_dict=wall_thickness_dict,
                                assessment_factor=3,
                                calculate_mean=True,)

save_df_pickle(filename='example_mean_monte-carlo_loop_df', df= df_mean, foldername='monte-carlo_output')

print('Mean exceedences:', round(len(df_mean.loc[df_mean.dw_concs > dw_norm]) / len(df_mean)*100, 3), '%, total sims:', len(df_mean) )

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
# -------------------------------------------------
# Plotting of the input parameters
# -------------------------------------------------
# Creating a series of data of in range of 1-10000.
save_results_to = check_create_folders(folder_name='figures')

seg1 = Segment(name='seg1',
            material= 'PE40',
            length=25,
            inner_diameter=0.1,
            wall_thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_conditions(
    chemical_name='Benzeen', 
                    concentration_groundwater =1,
                    concentration_drinking_water=0.1,
                    temperature_groundwater=12, 
                    flow_rate=0.5)


lps = np.arange(start=0, stop=10000, step=1)
cdf_K = []
cdf_D = []
cdf_K_age = []
cdf_D_age = []

# calculate the cdf for the different parameters...
for lp in lps:

    Sr_K = 0.31397266
    log_Kpw_ref = NormalDist(mu=seg1.log_Kpw_ref, #1.6476099999999998 #
                             sigma=Sr_K).inv_cdf(p=random.random())

    Sr_D = 0.19572320 #Sr = standard error of regression
    log_Dp_ref = NormalDist(mu=seg1.log_Dp_ref, #-11.54717333172 #
                            sigma=Sr_D).inv_cdf(p=random.random())
    
    Sr_age_D = 0.17 #Eqn 21, KWR 2016.056
    f_Dage = NormalDist(mu=0, 
                        sigma=Sr_age_D).inv_cdf(p=random.random())
    
    Sr_age_K = 0.05 #Eqn 22, KWR 2016.056
    f_Kage = NormalDist(mu=0, 
                        sigma=Sr_age_K).inv_cdf(p=random.random())

    cdf_K.append(log_Kpw_ref )
    cdf_D.append(log_Dp_ref )
    cdf_K_age.append(f_Kage)
    cdf_D_age.append(f_Dage )

cdf_K.sort()
cdf_D.sort()
cdf_K_age.sort()
cdf_D_age.sort()

#Plotting the Results
# plt.plot(cdf_K, lps, color = 'red')
plt.plot(cdf_D, lps/len(lps), color = 'blue')
plt.xlabel('log Dp ref.')
plt.ylabel('Cumulatieve dichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/log_Dp_ref.png', dpi=300, bbox_inches='tight')

#%%
plt.plot(cdf_K, lps/len(lps), color = 'blue')
plt.xlabel('log Kpw ref.')
plt.ylabel('Cumulatieve dichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/log_Kpw_ref.png', dpi=300, bbox_inches='tight')

#%%
plt.plot(cdf_K_age, lps/len(lps), color = 'blue')
plt.xlabel('f_k_age')
plt.ylabel('Cumulatieve dichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/f_k_age.png', dpi=300, bbox_inches='tight')
#%%
plt.plot(cdf_D_age, lps/len(lps), color = 'blue')
plt.xlabel('f_D_age')
plt.ylabel('Cumulatieve dichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/f_D_age.png', dpi=300, bbox_inches='tight')

# %%
save_results_to = check_create_folders(folder_name='figures')

df_mean_f_Dtemps = df_mean.sort_values(by='f_Dtemps')
df_mean_f_Dtemps.reset_index(inplace = True, drop=True)

plt.plot(df_mean_f_Dtemps.f_Dtemps, df_mean_f_Dtemps.index/len(df_mean), color = 'blue')
plt.xlabel('f_Dtemp')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/f_Dtemp.png', dpi=300, bbox_inches='tight')

#%%
save_results_to = check_create_folders(folder_name='figures')

dfx = df_mean.sort_values(by='activattion_energies')
dfx.reset_index(inplace = True, drop=True)

plt.plot(dfx.activattion_energies, dfx.index/len(df_mean), color = 'blue')
plt.xlabel('Activattion_energy')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/activattion_energies.png', dpi=300, bbox_inches='tight')

#%%
save_results_to = check_create_folders(folder_name='figures')

dfx = df_mean.sort_values(by='partitioning_enthalpies')
dfx.reset_index(inplace = True, drop=True)

plt.plot(dfx.partitioning_enthalpies, dfx.index/len(df_mean), color = 'blue')
plt.xlabel('Partitioning Enthalpie')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/partitioning_enthalpies.png', dpi=300, bbox_inches='tight')

#%%

save_results_to = check_create_folders(folder_name='figures')

dfx = df_mean.sort_values(by='f_Ktemps')
dfx.reset_index(inplace = True, drop=True)

plt.plot(dfx.f_Ktemps, dfx.index/len(df_mean), color = 'blue')
plt.xlabel('f_Ktemp')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/f_Ktemp.png', dpi=300, bbox_inches='tight')

#%%

save_results_to = check_create_folders(folder_name='figures')

dfx = df_mean.sort_values(by='gw_concs')
dfx.reset_index(inplace = True, drop=True)

plt.plot(dfx.gw_concs, dfx.index/len(df_mean), color = 'blue')
plt.xlabel('gw_concs')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/gw_concs.png', dpi=300, bbox_inches='tight')

#%%

save_results_to = check_create_folders(folder_name='figures')

dfx = df_mean.sort_values(by='f_Kconcs')
dfx.reset_index(inplace = True, drop=True)

plt.plot(dfx.f_Kconcs, dfx.index/len(df_mean), color = 'blue')
plt.xlabel('f_Kconcs')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/f_Kconcs.png', dpi=300, bbox_inches='tight')

#%%
save_results_to = check_create_folders(folder_name='figures')

dfx = df_mean.sort_values(by='f_Dconcs')
dfx.reset_index(inplace = True, drop=True)

plt.plot(dfx.f_Dconcs, dfx.index/len(df_mean), color = 'blue')
plt.xlabel('f_Dconcs')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/f_Dconcs.png', dpi=300, bbox_inches='tight')

#%%
df_plume_concs = load_pickle(filename='monte-carlo_plume_concs', foldername='monte-carlo_output')
df_plume_concs = df_plume_concs.sort_values(by='Extrapolated values')
df_plume_concs.reset_index(inplace = True, drop=True)

plt.plot(df_plume_concs['Extrapolated values'], df_plume_concs.index/len(df_plume_concs),)
plt.xlabel('Plume conc (mg/kg)')
plt.ylabel('Cumulatieve dichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/plume_conc.png', dpi=300, bbox_inches='tight')
