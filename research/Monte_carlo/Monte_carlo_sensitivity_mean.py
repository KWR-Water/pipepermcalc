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
def calculate_dw_concentration(parameter_dict, 
                               calculate_peak=False,
                                calculate_mean=False, 
                                column_name = None,
                                ):
    ''' 
    Calculate the peak or the mean dw concentration given the parameters in 
    the parameter_dict
    '''
    concentration_soil = parameter_dict['concentration_soil']
    pipe_length = parameter_dict['pipe_length']
    length_fraction_middle_point = parameter_dict['length_fraction_middle_point']
    length_plume = parameter_dict['length_plume']
    inner_diameter = parameter_dict['inner_diameter']
    flow_rate = parameter_dict['flow_rate']

    #calculate contact length of pipe w/contamination plume
    if length_fraction_middle_point >= length_plume / 2:
        contact_length = length_plume
    elif length_fraction_middle_point < length_plume / 2:
        contact_length = length_fraction_middle_point + length_plume / 2

    # @MartinvdS, do we vary this one too? or leave as the thickness of the median value?
    # AH_todo
    wall_thickness= (wall_thickness_dict[input_parameters['median_values']['inner_diameter']])

    assessment_factor = 3

    # Create pipe and set conditions
    # ------------------------------
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=contact_length,
                    inner_diameter=inner_diameter,
                    wall_thickness=wall_thickness)

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(concentration_soil = concentration_soil,
                        chemical_name="Benzeen", 
                        temperature_groundwater=12,
                        flow_rate=flow_rate, 
                        suppress_print=True)

    # Update the partitioning and diffusion coefficients
    # --------------------------------------------------

    seg1.log_Kpw_ref = parameter_dict['log_Kpw_ref']

    seg1.log_Dp_ref = parameter_dict['log_Dp_ref']

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

    #set assessment_factor
    pipe1.ASSESSMENT_FACTOR_GROUNDWATER = assessment_factor

    pipe1.validate_input_parameters()

    # Calculate concentrations, can we do in one loop and store seperate peak/mean conc
    if calculate_peak:
        dw_conc = pipe1.calculate_peak_dw_concentration()
    elif calculate_mean:
        dw_conc = pipe1.calculate_mean_dw_concentration()

    # Save data as df
    # ------------------------------
    data = {'concentration_drinking_water': pipe1.concentration_drinking_water, 
            'concentration_soil': pipe1.concentration_soil, 
            'concentration_groundwater': pipe1.concentration_groundwater, 
            'log_Kpw': seg1.log_Kpw, 
            'log_Dp': seg1.log_Dp, 
            'contact_length': seg1.length, 
            'pipe_length': pipe_length, 
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
            'partitioning_enthalpie': seg1.partitioning_enthalpie}

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

        df_data= calculate_dw_concentration(parameter_dict =  parameter_dict, 
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

# AH _todo -> convert these to simple pickle files 
# import the data 
# plume concentration
df = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')

# range lenth plume
plume_length_range = load_pickle(filename='monte_carlo_plume_lengths', foldername='inputs')
plume_length_values = list(plume_length_range['Diameter (m)'])

# location middle point: continuous equal distribution 0-1
length_fraction_middle_points = list(np.arange(start=0, stop=1.01, step=0.01))

# import data on pipe dimension and concentrations
df_plume_concs = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')
plume_concs = list(df_plume_concs['Extrapolated values'])

# Import the pipe information - inner diameter, thickness, length
df_PE40 = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40')

# range lenth pipe
df_PE40_Lengte_GIS = df_PE40.sort_values(by='Lengte_GIS')
pipe_length_values = list(df_PE40_Lengte_GIS.Lengte_GIS)

# range inner diameter
df_PE40_inner_diam = df_PE40[df_PE40['Binnendiam'] > 0]  
df_PE40_inner_diam = df_PE40_inner_diam.sort_values(by='Binnendiam')
df_PE40_inner_diam['Binnendiam'] = df_PE40_inner_diam['Binnendiam']  / 1000
inner_diam_values = list(df_PE40_inner_diam.Binnendiam)

# dictionary relating specific inner diameters to specific wall thicknesses    
wall_thickness_dict = {0.01239999962: 0.0018,
                        0.01560000038: 0.0022,
                        0.01960000038: 0.0027,
                        0.025: 0.0035,
                        0.03139999962: 0.0043,
                        0.03920000076: 0.0054,
                        0.04940000153: 0.0068, }

# AH_todo, finish this after disc. w/@Martin
# Water use per household => 1/3 = 1 person, 1/3 = 2 people, 1/9 each 3,4,5 people @ 125 L/p/d
# water_use = [0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.375, 0.500, 0.625] 
household_distribution = [1,1,1,1,1,1,1,1,1,1,
                            2,2,2,2,2,2,2,
                            3,3,3,
                            4,4,4,
                            5]

household_water_use = 0.125 # L/p/d
water_use =  (np.array(household_distribution) * household_water_use)

#%% Create dictionary with the median, st. devation of the input parameters
# -------------------------------------------------------------------------
input_parameters = {
            'median_values':
                    {'concentration_soil': df_plume_concs['Extrapolated values'].median(), 
                    'pipe_length':  df_PE40_Lengte_GIS.Lengte_GIS.median(),
                    'length_fraction_middle_point': float(np.median(length_fraction_middle_points)),
                    'length_plume':  plume_length_range['Diameter (m)'].median(),
                    'inner_diameter': df_PE40_inner_diam.Binnendiam.median(), 
                    'flow_rate': float(np.median(water_use)), 
                    'log_Dp_ref': -11.54717333172, 
                    'log_Kpw_ref': 1.6476099999999998, 
                    'DIFFUSION_A_C': 0.784077209735583, 
                    'PARTITIONING_A_C':  0.103965019849463,
                    'activattion_energy': 38.156061538172395, # for T = 12 deg. C
                    'partitioning_enthalpie': 8.943052711652054},# for T = 12 deg. C
                    
                'st_dev': 
                    {'concentration_soil': df_plume_concs['Extrapolated values'].std(), 
                    'pipe_length':  df_PE40_Lengte_GIS.Lengte_GIS.std(),
                    'length_fraction_middle_point': float(np.std(length_fraction_middle_points)),
                    'length_plume':  plume_length_range['Diameter (m)'].std(),
                    'inner_diameter': df_PE40_inner_diam.Binnendiam.std(),
                    'flow_rate': float(np.std(water_use)), 
                    'log_Dp_ref': 0.19572320, 
                    'log_Kpw_ref': 0.31397266, 
                    'DIFFUSION_A_C': 0.07662645, 
                    'PARTITIONING_A_C': 0.10106212,
                    'activattion_energy': 11.7958431, 
                    'partitioning_enthalpie':  13.2239059}}

# add dictionary with median + 1%
input_parameters['median_+1%'] = {}
for k,v in input_parameters['median_values'].items(): 
    input_parameters['median_+1%'][k] = float(v + np.abs(v *0.01))

    #@MartinvdS, check about 1% of LogDp_ref or 1% Dp_ref?? 
    # Makes problems for the st. dev calculations when changing so be careful

# add dictionary with median + 1 st. dev
input_parameters['median_+std'] = {}
for k,v in input_parameters['median_values'].items(): 
    input_parameters['median_+std'][k] = v  + input_parameters['st_dev'][k]

# override the values for the kow values
input_parameters['median_+1%'][ 'log_Dp_ref'] = np.log10(10**(input_parameters['median_values'][ 'log_Dp_ref']) + 10**(input_parameters['median_values'][ 'log_Dp_ref'])*0.01)
input_parameters['median_+1%'][ 'log_Kpw_ref'] = np.log10(10**(input_parameters['median_values'][ 'log_Kpw_ref']) + 10**(input_parameters['median_values'][ 'log_Kpw_ref'])*0.01)

#%%
save_results_to = check_create_folders(folder_name='output')

# Median values
# -------------
df_peak = calculate_dw_concentration(parameter_dict =  input_parameters['median_values'], 
                                     calculate_peak = True, 
                                     column_name = 'Median_peak')
df_mean = calculate_dw_concentration(parameter_dict =  input_parameters['median_values'], 
                                     calculate_mean = True,
                                     column_name = 'Median_mean')

df_peak.to_excel(save_results_to+"/peak_median.xlsx")
df_mean.to_excel(save_results_to+"/mean_median.xlsx")

#%%
# Median +1%
# ----------
df_peak_1 = calculate_dw_varying_one_parameter (option = 'median_+1%',
                                        input_parameters=input_parameters,
                                        calculate_peak=True,)
df_peak_1.to_excel(save_results_to+"/peak_1%.xlsx")
df_peak_1=df_peak_1.append(df_peak)
df_peak_1.insert(loc=0, column='C/Cm', value = df_peak_1['concentration_drinking_water'] / float(df_peak['concentration_drinking_water']))

df_mean_1 = calculate_dw_varying_one_parameter (option = 'median_+1%',
                                        input_parameters=input_parameters,
                                        calculate_mean=True, )
df_mean_1.to_excel(save_results_to+"/mean_1%.xlsx")
df_mean_1=df_mean_1.append(df_mean)
df_mean_1.insert(loc=0, column='C/Cm', value = df_mean_1['concentration_drinking_water'] / float(df_mean['concentration_drinking_water']))
#%%
# Median + 1 st. dev
# ------------------

df_peak_std = calculate_dw_varying_one_parameter (option = 'median_+std',
                                        input_parameters=input_parameters,
                                        calculate_peak=True,)
df_peak_std.to_excel(save_results_to+"/peak_1_st_dev.xlsx")
df_peak_std=df_peak_std.append(df_peak)
df_peak_std.insert(loc=0, column='C/Cm', value = df_peak_std['concentration_drinking_water'] / float(df_peak['concentration_drinking_water']))

df_mean_std = calculate_dw_varying_one_parameter (option = 'median_+std',
                                        input_parameters=input_parameters,
                                        calculate_mean=True, )
df_mean_std.to_excel(save_results_to+"/mean_1_st_dev.xlsx")
df_mean_std=df_mean_std.append(df_mean)
df_mean_std.insert(loc=0, column='C/Cm', value = df_mean_std['concentration_drinking_water'] / float(df_mean['concentration_drinking_water']))

# %%

df_inputs = pd.DataFrame(input_parameters)

# Mean
df_mean_summary = pd.DataFrame()

df_mean_summary['median + stdev']=df_mean_std['concentration_drinking_water']
df_mean_summary['median + stdev (%)']=df_mean_std['C/Cm']
df_mean_summary['median + 1%']=df_mean_1['concentration_drinking_water']
df_mean_summary['median + 1% (%)']=df_mean_1['C/Cm']

# Peak
df_peak_summary = pd.DataFrame()

df_peak_summary['median + stdev']=df_peak_std['concentration_drinking_water']
df_peak_summary['median + stdev (%)']=df_peak_std['C/Cm']
df_peak_summary['median + 1%']=df_peak_1['concentration_drinking_water']
df_peak_summary['median + 1% (%)']=df_peak_1['C/Cm']
df_peak_summary

#%%
# highlight values based on columns
df_peak_1=df_peak_1.style.background_gradient(axis=0)  
df_mean_1=df_mean_1.style.background_gradient(axis=0)  
df_peak_std=df_peak_std.style.background_gradient(axis=0)  
df_mean_std=df_mean_std.style.background_gradient(axis=0)  

with pd.ExcelWriter('output/sensitivity_results.xlsx') as writer:  
    df_inputs.to_excel(writer, sheet_name='Input_values')
    df_mean_summary.to_excel(writer, sheet_name='Mean_summary')
    df_peak_summary.to_excel(writer, sheet_name='Peak_summary')
    df_peak_1.to_excel(writer, sheet_name='peak_1%')
    df_mean_1.to_excel(writer, sheet_name='mean_1%')
    df_peak_std.to_excel(writer, sheet_name='peak_std')
    df_mean_std.to_excel(writer, sheet_name='mean_std')    

# left off: need to make a comparison table or graphic to show the influence of the different params
# %%
