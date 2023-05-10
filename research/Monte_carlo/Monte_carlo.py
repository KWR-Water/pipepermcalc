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

def run_Monte_Carlo_simulation (plume_concs, 
                                inner_diam_values,
                                wall_thickness_dict,
                                assessment_factor,
                                water_use,
                                pipe_length_values,
                                length_fraction_middle_point,
                                plume_length_values,
                                calculate_peak=False, 
                                calculate_mean=False,
                                simulations_per_batch = 1000,
                                tolerance = 0.01,
                                update_partitioning_coefficients = True,
                                ):

    # Loop through the combinations, save the dw concentrations
    dw_concs = []
    soil_concs = []
    gw_concs = []
    log_Kpws = []
    log_Dps = []
    flow_rates = []
    inner_diameters = []
    wall_thicknesses =[]

    length_pipes = []
    length_middle_points =[]
    length_plumes = []
    contact_lengths = []

    log_Dp_refs =[]
    log_Kpw_refs =[]

    a_c_Ds = []
    a_c_Ks = []
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
    # simulations_per_batch = range(1000)
    # tolerance = 0.01

    while True:

        for lp in tqdm(range(simulations_per_batch)):
            # Input variables
            # ---------------
            # Soil Concentration, NBNL data
            soil_conc = random.choice(plume_concs)

            # Length, PWN data given in meters
            length_pipe = random.choice(pipe_length_values)

            # uniform distribution of location of middle point plume
            length_middle_point = random.choice(length_fraction_middle_point) * length_pipe

            #NBNL #ah_todo
            length_plume = random.choice(plume_length_values)

            #calculate contact length of pipe w/contamination plume
            contact_length = min(length_pipe, length_plume, (length_plume / 2) + min ((length_pipe - length_middle_point), length_middle_point))

            # Inner diameter, PWN data, given in mm (convert to m)
            inner_diameter= random.choice(inner_diam_values) 

            # matching wall thicnkess from inner diameter, PWN data, given in mm (convert to m)
            wall_thickness = (wall_thickness_dict[inner_diameter])

            flow_rate = random.choice(water_use)

            # Create pipe and set conditions
            # ------------------------------
            seg1 = Segment(name='seg1',
                            material='PE40',
                            length=contact_length,
                            inner_diameter=inner_diameter,
                            wall_thickness=wall_thickness)

            pipe1 = Pipe(segment_list=[seg1])

            #set assessment_factor
            if type(assessment_factor) is int:
                pipe1.ASSESSMENT_FACTOR_GROUNDWATER = assessment_factor
            elif type(assessment_factor) is list:
                pipe1.ASSESSMENT_FACTOR_GROUNDWATER = random.choice(assessment_factor)

            pipe1.ASSESSMENT_FACTOR_SOIL = pipe1.ASSESSMENT_FACTOR_GROUNDWATER

            pipe1.set_conditions(concentration_soil = soil_conc,
                                chemical_name="Benzeen", 
                                temperature_groundwater=12,
                                flow_rate=flow_rate, 
                                suppress_print=True)

            if update_partitioning_coefficients:
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
                # Correction on the a_c factor
                Sr_conc_D = 0.07662645 #excel:'CONC' AE-4
                a_c_D =NormalDist(mu=seg1.DIFFUSION_A_C, #DIFFUSION_A_C = 0.784077209735583
                                    sigma=Sr_conc_D).inv_cdf(p=random.random()) 
                Cg_Sw = min((pipe1.concentration_groundwater/pipe1.solubility), 1)
                f_Dconc =  a_c_D * (Cg_Sw - seg1.DIFFUSION_CREF_SW) # DIFFUSION_CREF_SW = 0.5
                
                Sr_conc_K = 0.10106212 #excel:'CONC' W-4
                a_c_K =NormalDist(mu=seg1.PARTITIONING_A_C, #PARTITIONING_A_C = 0.103965019849463
                                    sigma=Sr_conc_K).inv_cdf(p=random.random()) 
                Cg_Sw = min((pipe1.concentration_groundwater/pipe1.solubility), 1)
                f_Kconc = a_c_K * (Cg_Sw - seg1.PARTITIONING_CREF_SW) # PARTITIONING_CREF_SW = 1.000
                
                # temperature corrections
                #corrections on the act. engergy/enthalpie itself not the factor
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
            
            else:
                log_Kpw_ref = seg1.log_Kpw_ref 
                a_c_K = seg1.PARTITIONING_A_C
                f_Kconc = seg1.f_Kconc
                activattion_energy = seg1.activattion_energy
                f_Ktemp = seg1.f_Ktemp
                f_Kage = 0


                log_Dp_ref = seg1.log_Dp_ref
                a_c_D = seg1.DIFFUSION_A_C
                f_Dconc = seg1.f_Dconc
                partitioning_enthalpie = seg1.partitioning_enthalpie
                f_Dtemp = seg1.f_Dtemp
                f_Dage = 0
            
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

            flow_rates.append(pipe1.flow_rate)
            inner_diameters.append(inner_diameter)
            wall_thicknesses.append(seg1.wall_thickness)

            length_pipes.append(length_pipe)
            length_middle_points.append(length_middle_point)
            length_plumes.append(length_plume)
            contact_lengths.append(seg1.length)

            log_Dp_refs.append(log_Dp_ref)
            log_Kpw_refs.append(log_Kpw_ref)

            a_c_Ds.append(a_c_D)
            a_c_Ks.append(a_c_K)
            f_Dconcs.append(f_Dconc)
            f_Kconcs.append(f_Kconc)

            activattion_energies.append(activattion_energy)
            partitioning_enthalpies.append(partitioning_enthalpie)
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
    data = zip(dw_concs, soil_concs, gw_concs, log_Kpws, log_Dps, 
            contact_lengths, length_pipes, length_middle_points, length_plumes,
            inner_diameters, flow_rates, wall_thicknesses, 
            log_Dp_refs, log_Kpw_refs, 
            f_Dconcs, f_Kconcs, a_c_Ds, a_c_Ks,
            f_Dtemps, f_Ktemps, activattion_energies, partitioning_enthalpies)
    df = pd.DataFrame(data,  
                      columns = ['dw_concs', 'soil_conc', 'gw_concs', 'Kpw', 'Dp', 
                    'contact_lengths', 'length_pipes', 'length_middle_points', 'length_plumes',
                    'inner_diameter', 'flow_rates', 'wall_thicknesses',
                    'log_Dp_refs', 'log_Kpw_refs',
                    'f_Dconcs', 'f_Kconcs', 'a_c_Ds', 'a_c_Ks',
                    'f_Dtemps', 'f_Ktemps', 'activattion_energies', 'partitioning_enthalpies'])
    df.sort_values(by=['dw_concs'], inplace=True)
    df.reset_index(inplace=True)

    return df
#%%
def plot_cumulative_distribution(df, dw_norm, save_name, save_results_to):
        fig = plt.figure(figsize=[10, 5])

        plt.plot(df.dw_concs, df.index/len(df), ) 
        plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
        plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
        plt.ylabel('Cumulative kansdichtheid')
        plt.title('Overschrijdingen per jaar: '+str(round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 1))+ '%, total sims:'+ str(len(df)) )
        plt.xscale('log')
        plt.xlim(1e-12, 1)
        plt.legend()
        plt.savefig(save_results_to +'/'+save_name + '.png', dpi=300, bbox_inches='tight')

#%%
def run_simulation_export_save_plot(dw_norm, 
                                input_parameters,
                                assessment_factor, 
                                calculate_mean, 
                                calculate_peak,
                                save_name, 
                                save_results_to, 
                                update_partitioning_coefficients = True,
                                simulations_per_batch = 1000 ):
    ''' 
    Run the Monte Carlo based on the input parameters dictionary, 
    export results to pickle and plot the exceedences

    '''
    
    df = run_Monte_Carlo_simulation (plume_concs=input_parameters['concentration_soil'], 
                                    inner_diam_values=input_parameters['inner_diameter'],
                                    wall_thickness_dict=input_parameters['wall_thickness_dict'],
                                    water_use=input_parameters['flow_rate'],
                                    pipe_length_values = input_parameters['length_pipe'],
                                    length_fraction_middle_point=input_parameters['length_fraction_middle_point'],
                                    plume_length_values=input_parameters['length_plume'],
                                    assessment_factor=assessment_factor,
                                    calculate_mean=calculate_mean,
                                    calculate_peak=calculate_peak,
                                    update_partitioning_coefficients=update_partitioning_coefficients,
                                    simulations_per_batch = simulations_per_batch)

    save_df_pickle(filename=save_name, df= df, foldername=save_results_to)
    df.to_excel(save_results_to+ "/" +save_name+".xlsx")

    plot_cumulative_distribution(df=df, 
                                dw_norm =dw_norm, 
                                save_name = save_name, 
                                save_results_to='figures')

    print('Mean exceedences:', round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 3), '%, total sims:', len(df) )
    return df
