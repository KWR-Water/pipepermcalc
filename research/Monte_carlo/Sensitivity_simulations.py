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

from Monte_carlo import *

#%%
def calculate_dw_concentration( wall_thickness_dict,
                                parameter_dict, 
                               calculate_peak=False,
                                calculate_mean=False, 
                                column_name = None,
                                assessment_factor = 3
                                ):
    ''' 
    Calculate the peak or the mean dw concentration given the parameters in 
    the parameter_dict
    '''
    concentration_soil = parameter_dict['concentration_soil']
    length_pipe = parameter_dict['length_pipe']
    length_fraction_middle_point = parameter_dict['length_fraction_middle_point']
    length_plume = parameter_dict['length_plume']
    inner_diameter = parameter_dict['inner_diameter']
    flow_rate = parameter_dict['flow_rate']

    length_middle_point = length_fraction_middle_point * length_pipe

    #calculate contact length of pipe w/contamination plume
    contact_length = min(length_pipe, length_plume, (length_plume / 2) + min ((length_pipe - length_middle_point), length_middle_point))

    wall_thickness= (wall_thickness_dict[input_parameters['median_values']['inner_diameter']])

    # Create pipe and set conditions
    # ------------------------------
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=contact_length,
                    inner_diameter=inner_diameter,
                    wall_thickness=wall_thickness)

    pipe1 = Pipe(segment_list=[seg1])

    #set assessment_factor
    pipe1.ASSESSMENT_FACTOR_GROUNDWATER = assessment_factor

    pipe1.ASSESSMENT_FACTOR_SOIL = pipe1.ASSESSMENT_FACTOR_GROUNDWATER

    pipe1.set_conditions(concentration_soil = concentration_soil,
                        chemical_name="Benzeen", 
                        temperature_groundwater=12,
                        flow_rate=flow_rate, 
                        suppress_print=True)

    # Update the partitioning and diffusion coefficients
    # --------------------------------------------------

    seg1.log_Kpw_ref = parameter_dict['log_Kpw_ref']

    seg1.log_Dp_ref = parameter_dict['log_Dp_ref']
    
    pipe1.log_distribution_coefficient = parameter_dict['log_distribution_coefficient']
    pipe1.concentration_groundwater = ((pipe1.concentration_soil * pipe1.ASSESSMENT_FACTOR_GROUNDWATER) 
                                                / ( 10 ** pipe1.log_distribution_coefficient * pipe1.ASSESSMENT_FACTOR_SOIL ))

    seg1.DIFFUSION_A_C = parameter_dict['DIFFUSION_A_C']
    Cg_Sw = min((pipe1.concentration_groundwater / pipe1.solubility), 1)
    seg1.f_Dconc = seg1.DIFFUSION_A_C * (Cg_Sw - seg1.DIFFUSION_CREF_SW) # DIFFUSION_CREF_SW = 0.5

    seg1.PARTITIONING_A_C = parameter_dict['PARTITIONING_A_C']
    Cg_Sw = min((pipe1.concentration_groundwater / pipe1.solubility), 1)
    seg1.f_Kconc = seg1.PARTITIONING_A_C  * (Cg_Sw - seg1.PARTITIONING_CREF_SW) # PARTITIONING_CREF_SW = 1.000

    seg1.activattion_energy = parameter_dict['activattion_energy'] 
    seg1.f_Dtemp = (seg1.activattion_energy / (0.008314 * np.log(10)) 
            * (1 / (25 + 273) - 1 / (pipe1.temperature_groundwater + 273)))
    
    seg1.partitioning_enthalpie = parameter_dict['partitioning_enthalpie']
    seg1.f_Ktemp = ( seg1.partitioning_enthalpie / (0.008314 * np.log(10)) 
            * (1 / (25 + 273) - 1 / (pipe1.temperature_groundwater + 273)))

    f_Kage = 0
    f_Dage = 0

    # Set the Kpw and Dp
    seg1.log_Kpw = seg1.log_Kpw_ref + seg1.f_Kconc + seg1.f_Ktemp + f_Kage
    seg1.log_Dp = seg1.log_Dp_ref + seg1.f_Dconc + seg1.f_Dtemp + f_Dage

    pipe1.validate_input_parameters()

    # Calculate concentrations, can we do in one loop and store seperate peak/mean conc
    if calculate_peak:
        seg1.stagnation_factor = parameter_dict['stagnation_factor']
        dw_conc = pipe1.calculate_peak_dw_concentration()
    elif calculate_mean:
        seg1.stagnation_factor = 0
        dw_conc = pipe1.calculate_mean_dw_concentration()
     

    # Save data as df
    # ------------------------------
    data = {'concentration_drinking_water': pipe1.concentration_drinking_water, 
            'concentration_soil': pipe1.concentration_soil, 
            'concentration_groundwater': pipe1.concentration_groundwater, 
            'log_Kpw': seg1.log_Kpw, 
            'log_Dp': seg1.log_Dp, 
            'contact_length': seg1.length, 
            'length_pipe': length_pipe, 
            'length_fraction_middle_point': length_fraction_middle_point, 
            'length_plume': length_plume,
            'inner_diameter': seg1.inner_diameter, 
            'flow_rate': pipe1.flow_rate, 
            'wall_thickness': seg1.wall_thickness, 
            'log_Dp_ref': seg1.log_Dp_ref, 
            'log_Kpw_ref': seg1.log_Kpw_ref, 
            'f_Dconc': seg1.f_Dconc, 
            'f_Kconc': seg1.f_Kconc, 
            'DIFFUSION_A_C': seg1.DIFFUSION_A_C, 
            'PARTITIONING_A_C': seg1.PARTITIONING_A_C,
            'f_Dtemp': seg1.f_Dtemp, 
            'f_Ktemp': seg1.f_Ktemp, 
            'activattion_energy': seg1.activattion_energy, 
            'partitioning_enthalpie': seg1.partitioning_enthalpie, 
            'log_distribution_coefficient': pipe1.log_distribution_coefficient, 
            'stagnation_factor': seg1.stagnation_factor
            }

    df = pd.DataFrame(data, index=[column_name])

    return df

def calculate_dw_varying_one_parameter (option, 
                                        input_parameters, 
                                        calculate_peak=False,
                                        calculate_mean=False, 
                                        ):
    variable_parameters = list(input_parameters['median_values'].keys())
    dfs = []

    for variable_parameter in variable_parameters:
        parameter_dict = {}

        for var in variable_parameters:
            if var == variable_parameter:
                parameter_dict[var] = input_parameters[option][var]
            else:
                parameter_dict[var] = input_parameters['median_values'][var]

        df_data= calculate_dw_concentration(wall_thickness_dict = input_parameters['cumulative_distributions']['wall_thickness_dict'], 
                                            parameter_dict =  parameter_dict, 
                                            calculate_peak=calculate_peak,
                                            calculate_mean=calculate_mean, 
                                            # column_name=option+'_'+variable_parameter)
                                            column_name=variable_parameter)
        dfs.append(df_data)

    df = pd.concat(dfs)
    return df

#%%
# Input Parameters
# ----------------
def create_input_parameters_dict ():
    ''' 
    Create input parameters dictionary with the cumulative distributions, 
    median, standard deviation and % variations
    '''

    input_parameters = {}
    input_parameters['cumulative_distributions'] = load_pickle(filename='monte_carlo_input_parameters', foldername='inputs')

    # Create dictionary with the median, st. devation of the input parameters
    median_values = {'concentration_soil': float(np.median(input_parameters['cumulative_distributions'] ['concentration_soil'])), 
                        'length_pipe': float(np.median(input_parameters['cumulative_distributions'] ['length_pipe'])),
                        'length_fraction_middle_point': float(np.median(input_parameters['cumulative_distributions'] ['length_fraction_middle_point'])),
                        'length_plume':  float(np.median(input_parameters['cumulative_distributions'] ['length_plume'])),
                        'inner_diameter': float(np.median(input_parameters['cumulative_distributions'] ['inner_diameter'])), 
                        'flow_rate': float(np.median(input_parameters['cumulative_distributions'] ['flow_rate'])), 
                        'log_Dp_ref': -11.54717333172, 
                        'log_Kpw_ref': 1.6476099999999998, 
                        'DIFFUSION_A_C': 0.784077209735583, 
                        'PARTITIONING_A_C':  0.103965019849463,
                        'activattion_energy': 38.156061538172395, # for T = 12 deg. C
                        'partitioning_enthalpie': 8.943052711652054, # for T = 12 deg. C
                        'log_distribution_coefficient': 0.659555885,
                        'stagnation_factor': 1.3868499849384468
                        }
                        
    st_dev = {'concentration_soil': float(np.std(input_parameters['cumulative_distributions'] ['concentration_soil'])), 
                        'length_pipe': float(np.std(input_parameters['cumulative_distributions'] ['length_pipe'])),
                        'length_fraction_middle_point': float(np.std(input_parameters['cumulative_distributions'] ['length_fraction_middle_point'])),
                        'length_plume':  float(np.std(input_parameters['cumulative_distributions'] ['length_plume'])),
                        'inner_diameter': float(np.std(input_parameters['cumulative_distributions'] ['inner_diameter'])), 
                        'flow_rate': float(np.std(input_parameters['cumulative_distributions'] ['flow_rate'])),
                        'log_Dp_ref': 0.19572320, 
                        'log_Kpw_ref': 0.31397266, 
                        'DIFFUSION_A_C': 0.07662645, 
                        'PARTITIONING_A_C': 0.10106212,
                        'activattion_energy': 11.7958431, 
                        'partitioning_enthalpie':  13.2239059,
                        'log_distribution_coefficient': 0,
                        'stagnation_factor': 0
                        }

    input_parameters['median_values'] = median_values
    input_parameters['st_dev'] = st_dev

    percent_change = 0.1

    # add dictionary with median + %
    input_parameters['median_+%'] = {}
    for k,v in input_parameters['median_values'].items(): 
        input_parameters['median_+%'][k] = float(v + np.abs(v * percent_change))

        #@MartinvdS, check about % of LogDp_ref or % Dp_ref?? 
        # Makes problems for the st. dev calculations when changing so be careful

    # add dictionary with median + 1 st. dev
    input_parameters['median_+std'] = {}
    for k,v in input_parameters['median_values'].items(): 
        input_parameters['median_+std'][k] = v  + input_parameters['st_dev'][k]

    # override the values for the kpw, Dp, and Kbw values
    input_parameters['median_+%'][ 'log_Dp_ref'] = np.log10(10**(input_parameters['median_values'][ 'log_Dp_ref']) + 10**(input_parameters['median_values'][ 'log_Dp_ref']) * percent_change)
    input_parameters['median_+%'][ 'log_Kpw_ref'] = np.log10(10**(input_parameters['median_values'][ 'log_Kpw_ref']) + 10**(input_parameters['median_values'][ 'log_Kpw_ref']) * percent_change)
    input_parameters['median_+%'][ 'log_distribution_coefficient'] = 0.700948570158225 #np.log10(10**(input_parameters['median_values'][ 'log_distribution_coefficient']) + 10**(input_parameters['median_values'][ 'log_distribution_coefficient']) * percent_change)

    return input_parameters

input_parameters = create_input_parameters_dict ()
#%%
save_results_to = check_create_folders(folder_name='output')

# Median values
# -------------
df_peak = calculate_dw_concentration(wall_thickness_dict = input_parameters['cumulative_distributions']['wall_thickness_dict'], 
                                    parameter_dict =  input_parameters['median_values'], 
                                     calculate_peak = True, 
                                     column_name = 'Median_peak')
df_mean = calculate_dw_concentration(wall_thickness_dict = input_parameters['cumulative_distributions']['wall_thickness_dict'], 
                                     parameter_dict =  input_parameters['median_values'], 
                                     calculate_mean = True,
                                     column_name = 'Median_mean')
df_peak['concentration_drinking_water'] *1000*1000, df_mean['concentration_drinking_water'] *1000*1000

#%%
# Median values, assessment_factor = 1
# -------------
df_peak_f1 = calculate_dw_concentration(wall_thickness_dict = input_parameters['cumulative_distributions']['wall_thickness_dict'], 
                                    parameter_dict =  input_parameters['median_values'], 
                                     calculate_peak = True, 
                                     column_name = 'Median_peak', 
                                     assessment_factor=1)

df_mean_f1 = calculate_dw_concentration(wall_thickness_dict = input_parameters['cumulative_distributions']['wall_thickness_dict'], 
                                     parameter_dict =  input_parameters['median_values'], 
                                     calculate_mean = True,
                                     column_name = 'Median_mean', 
                                      assessment_factor=1)

df_peak_f1['concentration_drinking_water'] *1000*1000, df_mean_f1['concentration_drinking_water'] *1000*1000
#%%
# Median +%
# ----------
df_peak_1 = calculate_dw_varying_one_parameter (option = 'median_+%',
                                        input_parameters=input_parameters,
                                        calculate_peak=True,)
# df_peak_1.to_excel(save_results_to+"/peak_%.xlsx")
df_peak_1=df_peak_1.append(df_peak)
df_peak_1.insert(loc=0, column='C/Cm', value = df_peak_1['concentration_drinking_water'] / float(df_peak['concentration_drinking_water']))

df_mean_1 = calculate_dw_varying_one_parameter (option = 'median_+%',
                                        input_parameters=input_parameters,
                                        calculate_mean=True, )
# df_mean_1.to_excel(save_results_to+"/mean_%.xlsx")
df_mean_1=df_mean_1.append(df_mean)
df_mean_1.insert(loc=0, column='C/Cm', value = df_mean_1['concentration_drinking_water'] / float(df_mean['concentration_drinking_water']))
#%%
# Median + 1 st. dev
# ------------------

df_peak_std = calculate_dw_varying_one_parameter (option = 'median_+std',
                                        input_parameters=input_parameters,
                                        calculate_peak=True,)
# df_peak_std.to_excel(save_results_to+"/peak_1_st_dev.xlsx")
df_peak_std=df_peak_std.append(df_peak)
df_peak_std.insert(loc=0, column='C/Cm', value = df_peak_std['concentration_drinking_water'] / float(df_peak['concentration_drinking_water']))

df_mean_std = calculate_dw_varying_one_parameter (option = 'median_+std',
                                        input_parameters=input_parameters,
                                        calculate_mean=True, )
# df_mean_std.to_excel(save_results_to+"/mean_1_st_dev.xlsx")
df_mean_std=df_mean_std.append(df_mean)
df_mean_std.insert(loc=0, column='C/Cm', value = df_mean_std['concentration_drinking_water'] / float(df_mean['concentration_drinking_water']))

# %%
# Create summary files for the inputs, mean and peak simulations
input_parameters.pop('cumulative_distributions')
df_inputs = pd.DataFrame(input_parameters)
df_inputs['peak_dw_conc'] = float(df_peak['concentration_drinking_water'])
df_inputs['mean_dw_conc'] = float(df_mean['concentration_drinking_water'])

# Mean
df_mean_summary = pd.DataFrame()

df_mean_summary['median + stdev']=df_mean_std['concentration_drinking_water']
df_mean_summary['median + stdev (%)']=df_mean_std['C/Cm']
df_mean_summary['median + %']=df_mean_1['concentration_drinking_water']
df_mean_summary['median + % (%)']=df_mean_1['C/Cm']

# Peak
df_peak_summary = pd.DataFrame()

df_peak_summary['median + stdev']=df_peak_std['concentration_drinking_water']
df_peak_summary['median + stdev (%)']=df_peak_std['C/Cm']
df_peak_summary['median + %']=df_peak_1['concentration_drinking_water']
df_peak_summary['median + % (%)']=df_peak_1['C/Cm']
df_peak_summary

df_analysis = pd.DataFrame()
df_analysis['(median+std)/median'] = df_inputs['median_+std'] / df_inputs['median_values']
df_analysis['mean_median + %'] = df_mean_summary['median + % (%)']
df_analysis['mean_median + stdev']= df_mean_summary['median + stdev (%)']
df_analysis['peak_median + %'] = df_peak_summary['median + % (%)']
df_analysis['peak_median + stdev']= df_peak_summary['median + stdev (%)']

# override the values for the Kpw and Dp since they don't calculate the same way
df_analysis['(median+std)/median'].loc['log_Dp_ref'] =  10**(df_inputs['median_+std'].loc['log_Dp_ref']) / 10**(df_inputs['median_values'].loc['log_Dp_ref'])
df_analysis['(median+std)/median'].loc['log_Kpw_ref'] = 10**(df_inputs['median_+std'].loc['log_Kpw_ref']) / 10**(df_inputs['median_values'].loc['log_Kpw_ref'])
# df_analysis['(median+std)/median'].loc['log_distribution_coefficient'] = 10**(df_inputs['median_+std'].loc['log_distribution_coefficient']) / 10**(df_inputs['median_values'].loc['log_distribution_coefficient'])

df_analysis['mean_median + %_norm'] = np.abs(df_mean_summary['median + % (%)'] - 1)
df_analysis['mean_median + stdev_norm'] = np.abs(df_mean_summary['median + stdev (%)'] - 1)
df_analysis['peak_median + %_norm'] = np.abs(df_peak_summary['median + % (%)'] - 1)
df_analysis['peak_median + stdev_norm'] = np.abs(df_peak_summary['median + stdev (%)'] - 1)

#%%
# highlight values based on columns
df_peak_1=df_peak_1.style.background_gradient(axis=0)  
df_mean_1=df_mean_1.style.background_gradient(axis=0)  
df_peak_std=df_peak_std.style.background_gradient(axis=0)  
df_mean_std=df_mean_std.style.background_gradient(axis=0)  

with pd.ExcelWriter('output/sensitivity_results.xlsx') as writer:  
    df_analysis.to_excel(writer, sheet_name='df_analysis')
    df_inputs.to_excel(writer, sheet_name='Input_values')
    df_mean_summary.to_excel(writer, sheet_name='Mean_summary')
    df_peak_summary.to_excel(writer, sheet_name='Peak_summary')
    df_peak_1.to_excel(writer, sheet_name='peak_%')
    df_mean_1.to_excel(writer, sheet_name='mean_%')
    df_peak_std.to_excel(writer, sheet_name='peak_std')
    df_mean_std.to_excel(writer, sheet_name='mean_std')    

# left off: need to make a comparison table or graphic to show the influence of the different params
# %%
