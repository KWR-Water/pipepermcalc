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
# Import the input parameter data
input_parameters = load_pickle(filename='monte_carlo_input_parameters', foldername='inputs')

#%%
save_results_to = check_create_folders(folder_name='monte_carlo_output')

run_sim = True
if run_sim:
    # Assessment factor = 3
    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 3, 
                        calculate_mean = True, 
                        calculate_peak=False,
                        save_name = 'mean_monte_carlo_F_as=3', 
                        save_results_to=save_results_to)

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 3, 
                        calculate_mean = False, 
                        calculate_peak=True,
                        save_name = 'peak_monte_carlo_F_as=3', 
                        save_results_to=save_results_to)


    # Assessment factor = 1
    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 1, 
                        calculate_mean = True, 
                        calculate_peak=False,
                        save_name = 'mean_monte_carlo_F_as=1', 
                        save_results_to=save_results_to)

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 1, 
                        calculate_mean = False, 
                        calculate_peak=True,
                        save_name = 'peak_monte_carlo_F_as=1', 
                        save_results_to=save_results_to)

    # Assessment factor = 7
    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 7, 
                        calculate_mean = True, 
                        calculate_peak=False,
                        save_name = 'mean_monte_carlo_F_as=7', 
                        save_results_to=save_results_to)

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 7, 
                        calculate_mean = False, 
                        calculate_peak=True,
                        save_name = 'peak_monte_carlo_F_as=7', 
                            save_results_to=save_results_to)
    # Fass = 3
    # Csoil & Lplume = cumulative distribution
    # Other parameters = median
    input_parameters = load_pickle(filename='monte_carlo_input_parameters', foldername='inputs')

    input_parameters.update({'length_pipe': [float(np.median(input_parameters['length_pipe'])), float(np.median(input_parameters['length_pipe']))],
                        'length_fraction_middle_point': [float(np.median(input_parameters ['length_fraction_middle_point'])), float(np.median(input_parameters ['length_fraction_middle_point']))],
                        'inner_diameter': [float(np.median(input_parameters['inner_diameter'])), float(np.median(input_parameters['inner_diameter']))], 
                        'flow_rate': [float(np.median(input_parameters ['flow_rate'])), float(np.median(input_parameters ['flow_rate']))] })

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 3, 
                        calculate_mean = True, 
                        calculate_peak=False,
                        save_name = 'mean_monte_carlo_F_as=3_median', 
                        save_results_to=save_results_to, 
                        update_partitioning_coefficients = False, )

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 3, 
                        calculate_mean = False, 
                        calculate_peak=True,
                        save_name = 'peak_monte_carlo_F_as=3_median', 
                        save_results_to=save_results_to,
                        update_partitioning_coefficients = False,)

    # Fass = 1
    # Csoil & Lplume = cumulative distribution
    # Other parameters = median
    input_parameters = load_pickle(filename='monte_carlo_input_parameters', foldername='inputs')

    input_parameters.update({'length_pipe': [float(np.median(input_parameters['length_pipe'])), float(np.median(input_parameters['length_pipe']))],
                        'length_fraction_middle_point': [float(np.median(input_parameters ['length_fraction_middle_point'])), float(np.median(input_parameters ['length_fraction_middle_point']))],
                        'inner_diameter': [float(np.median(input_parameters['inner_diameter'])), float(np.median(input_parameters['inner_diameter']))], 
                        'flow_rate': [float(np.median(input_parameters ['flow_rate'])), float(np.median(input_parameters ['flow_rate']))] })

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 1, 
                        calculate_mean = True, 
                        calculate_peak=False,
                        save_name = 'mean_monte_carlo_F_as=1_median', 
                        save_results_to=save_results_to, 
                        update_partitioning_coefficients = False, )

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = 1, 
                        calculate_mean = False, 
                        calculate_peak=True,
                        save_name = 'peak_monte_carlo_F_as=1_median', 
                        save_results_to=save_results_to,
                        update_partitioning_coefficients = False,)

    
    # Fass = cumulative distribution 
    # Csoil & Lplume = cumulative distribution
    # Other parameters = median
    input_parameters = load_pickle(filename='monte_carlo_input_parameters', foldername='inputs')

    input_parameters.update({'length_pipe': [float(np.median(input_parameters['length_pipe'])), float(np.median(input_parameters['length_pipe']))],
                        'length_fraction_middle_point': [float(np.median(input_parameters ['length_fraction_middle_point'])), float(np.median(input_parameters ['length_fraction_middle_point']))],
                        'inner_diameter': [float(np.median(input_parameters['inner_diameter'])), float(np.median(input_parameters['inner_diameter']))], 
                        'flow_rate': [float(np.median(input_parameters ['flow_rate'])), float(np.median(input_parameters ['flow_rate']))] })

    # From '20190621 overschrijdingskans obv evaluatie praktijkmetingen.xlsx'
    assessment_factor = [4.670,7.162,4.071,76.460,409.059,2.293,4.029,61.222,209.175,3.129]

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = assessment_factor, 
                        calculate_mean = True, 
                        calculate_peak=False,
                        save_name = 'mean_monte_carlo_F_as=cum_dist_median', 
                        save_results_to=save_results_to, 
                        update_partitioning_coefficients = False, )

    df = run_simulation_export_save_plot(dw_norm = 0.001, 
                        input_parameters = input_parameters,
                        assessment_factor = assessment_factor, 
                        calculate_mean = False, 
                        calculate_peak=True,
                        save_name = 'peak_monte_carlo_F_as=cum_dist_median', 
                        save_results_to=save_results_to,
                        update_partitioning_coefficients = False,)
#%%

# save_results_to = check_create_folders(folder_name='monte_carlo_output')

# df_peak = load_pickle(filename='peak_monte_carlo_loop_df', foldername='monte_carlo_output')

# df_peak=df_peak.style.background_gradient(axis=0)  
# df_peak.to_excel(save_results_to+ "/peak_monte_carlo_loop_df.xlsx")

# df_mean = load_pickle(filename='mean_monte_carlo_loop_df', foldername='monte_carlo_output')
# df_mean=df_mean.style.background_gradient(axis=0)  
# df_mean.to_excel(save_results_to+ "/mean_monte_carlo_loop_df.xlsx")

#%%
import seaborn as sns

df_mean = load_pickle(filename='mean_monte_carlo_F_as=3', foldername='monte_carlo_output')
# calculate the correlation matrix
mean_corr = df_mean.corr()
mean_corr1 = mean_corr[['dw_concs']]
mean_corr1.drop(['index','dw_concs',  'gw_concs', 'Kpw', 'Dp'], axis = 'index', inplace = True)

# plot the heatmap
fig = plt.figure(figsize=[10, 5])
sns.heatmap(mean_corr1, 
        # xticklabels=corr.columns,
        yticklabels=mean_corr1.index, vmin=0, vmax=0.5)
plt.savefig('figures/mean_monte_carlo_correlation.png', dpi=300, bbox_inches='tight')

df_peak = load_pickle(filename='mean_monte_carlo_F_as=3', foldername='monte_carlo_output')

# calculate the correlation matrix
peak_corr = df_peak.corr()
peak_corr1 = mean_corr[['dw_concs']]
peak_corr1.drop(['index','dw_concs',  'gw_concs', 'Kpw', 'Dp'], axis = 'index', inplace = True)

# plot the heatmap
fig = plt.figure(figsize=[10, 5])
sns.heatmap(peak_corr1, 
        # xticklabels=corr.columns,
        yticklabels=peak_corr1.index, vmin=0, vmax=0.5)
plt.savefig('figures/peak_monte_carlo_correlation.png', dpi=300, bbox_inches='tight')

with pd.ExcelWriter('monte_carlo_output/mean_monte_carlo_F_as=3_results_and_correlation.xlsx') as writer:  
    df_peak.to_excel(writer, sheet_name='df_peak')
    peak_corr1.to_excel(writer, sheet_name='peak_correlation')
    df_mean.to_excel(writer, sheet_name='df_mean')
    mean_corr1.to_excel(writer, sheet_name='Mean_correlation')

#%%
# df_mean = load_pickle(filename='mean_monte_carlo_loop_df', foldername='monte_carlo_output')
df_plot = df_mean
fig = plt.figure(figsize=[10, 5])
sns.scatterplot(data = df_plot, y = 'soil_conc', x = 'dw_concs', hue = 'log_Kpw_refs', linewidth=0 )
plt.xscale('log')
plt.yscale('log')
plt.savefig('figures/mean_monte_carlo_soil_v_dw_hue_Kpw_ref.png', dpi=300, bbox_inches='tight')
#%%

fig = plt.figure(figsize=[10, 5])
sns.scatterplot(data = df_plot, y = 'soil_conc', x = 'dw_concs', hue = 'log_Dp_refs', linewidth=0 )
plt.xscale('log')
plt.yscale('log')
plt.savefig('figures/mean_monte_carlo_soil_v_dw_hue_Dp_ref.png', dpi=300, bbox_inches='tight')

#%%
fig = plt.figure(figsize=[10, 5])
sns.scatterplot(data = df_plot, y = 'log_Kpw_refs', x = 'dw_concs', linewidth=0, hue = 'soil_conc',)
# plt.yscale('log')
plt.xscale('log')

plt.savefig('figures/mean_monte_carlo_Kpw_ref_v_dw.png', dpi=300, bbox_inches='tight')

#%%
fig = plt.figure(figsize=[10, 5])
sns.scatterplot(data = df_plot, y = 'log_Dp_refs', x = 'dw_concs', linewidth=0, hue = 'soil_conc',)
# plt.yscale('log')
plt.xscale('log')
plt.savefig('figures/mean_monte_carlo_Dp_ref_v_dw.png', dpi=300, bbox_inches='tight')

# fig = plt.figure(figsize=[10, 5])
# sns.scatterplot(data = df_plot, x = 'soil_conc', y = 'dw_concs', hue = 'contact_lengths')
#%%