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
# import the data on plume concentration
# df = pd.read_excel('inputs/OS fase 1 voor KWR.xlsx', 
#                     sheet_name='ah_data', header=[3], usecols = "N", ) 
# df = df[:144]
# save_df_pickle(filename='monte_carlo_plume_concs', df= df, foldername='inputs')
df = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')

# # range lenth plume
# plume_length_range =pd.read_excel('inputs/NBNL_plume_SA.xlsx', 
#                     sheet_name='plume_SA', header=[3], usecols = "U", nrows=12) 
# save_df_pickle(filename='monte_carlo_plume_lengths', df= plume_length_range, foldername='inputs')
plume_length_range = load_pickle(filename='monte_carlo_plume_lengths', foldername='inputs')
plume_length_values = list(plume_length_range['Diameter (m)'])

# location middle point: continuous equal distribution 0-1
length_fraction_middle_point = list(np.arange(start=0, stop=1, step=0.01))

# import data on pipe dimension and concentrations
df = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')
plume_concs = list(df['Extrapolated values'])

# Import the pipe information - inner diameter, thickness, length
df_PE40 = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40')

# range lenth pipe
df_PE40_Lengte_GIS = df_PE40.sort_values(by='Lengte_GIS')
pipe_length_values = list(df_PE40_Lengte_GIS.Lengte_GIS)

# range inner diameter
df_PE40_inner_diam = df_PE40[df_PE40['Binnendiam'] > 0]  
df_PE40_inner_diam = df_PE40_inner_diam.sort_values(by='Binnendiam')
df_PE40_inner_diam['Binnendiam'] = df_PE40_inner_diam['Binnendiam']  / 1000
df_PE40_inner_diam = df_PE40_inner_diam.round({'Binnendiam': 4})
inner_diam_values = list(df_PE40_inner_diam.Binnendiam)

# dictionary relating specific inner diameters to specific wall thicknesses    
wall_thickness_dict = {0.0124: 0.0018,
                        0.0156: 0.0022,
                        0.0196: 0.0027,
                        0.0250: 0.0035,
                        0.0314: 0.0043,
                        0.0392: 0.0054,
                        0.0494: 0.0068, }

# AH_todo, finish this after disc. w/@Martin
# Water use per household => 1/3 = 1 person, 1/3 = 2 people, 1/9 each 3,4,5 people @ 125 L/p/d
# water_use = [0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.375, 0.500, 0.625] 
household_distribution = [1,1,1,1,1,1,1,1,1,1,
                            2,2,2,2,2,2,2,
                            3,3,3,
                            4,4,4,
                            5]

household_water_use = 0.125 # L/p/d
water_use =  (np.array(household_distribution) * household_water_use).tolist()
#%%
save_results_to = check_create_folders(folder_name='monte_carlo_output')

df_peak = run_Monte_Carlo_simulation (plume_concs=plume_concs, 
                                inner_diam_values=inner_diam_values,
                                wall_thickness_dict=wall_thickness_dict,
                                water_use=water_use,
                                pipe_length_values = pipe_length_values,
                                length_fraction_middle_point=length_fraction_middle_point,
                                plume_length_values=plume_length_values,
                                assessment_factor=3,
                                calculate_peak=True,)

save_df_pickle(filename='peak_monte_carlo_loop_df', df= df_peak, foldername='monte_carlo_output')
df_peak.to_excel(save_results_to+ "/peak_monte_carlo_loop_df.xlsx")

dw_norm = 0.001 # Benzene drinking water norm, g/m3
print('Peak exceedences:', round(len(df_peak.loc[df_peak.dw_concs > dw_norm]) / len(df_peak)*100, 3), '%, total sims:', len(df_peak) )

df_mean = run_Monte_Carlo_simulation (plume_concs=plume_concs, 
                                inner_diam_values=inner_diam_values,
                                wall_thickness_dict=wall_thickness_dict,
                                water_use=water_use,
                                pipe_length_values = pipe_length_values,
                                length_fraction_middle_point=length_fraction_middle_point,
                                plume_length_values=plume_length_values,
                                assessment_factor=3,
                                calculate_mean=True,)

save_df_pickle(filename='mean_monte_carlo_loop_df', df= df_mean, foldername='monte_carlo_output')
df_mean.to_excel(save_results_to+ "/mean_monte_carlo_loop_df.xlsx")

print('Mean exceedences:', round(len(df_mean.loc[df_mean.dw_concs > dw_norm]) / len(df_mean)*100, 3), '%, total sims:', len(df_mean) )
#%%
df_mean.describe()

#%%
dw_norm = 0.001
save_results_to = check_create_folders(folder_name='figures')

df = load_pickle(filename='mean_monte_carlo_loop_df', foldername='monte_carlo_output')
fig = plt.figure(figsize=[10, 5])

plt.plot(df.dw_concs, df.index/len(df), ) 
plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
plt.ylabel('Cumulative kansdichtheid')
plt.title('Overschrijdingen per jaar: '+str(round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 1))+ '%, total sims:'+ str(len(df)) )
plt.xscale('log')
plt.xlim(1e-12, 1)
plt.legend()
plt.savefig(save_results_to +'/mean_monte_carlo_simulation.png', dpi=300, bbox_inches='tight')

#%%
save_results_to = check_create_folders(folder_name='figures')

df = load_pickle(filename='peak_monte_carlo_loop_df', foldername='monte_carlo_output')
fig = plt.figure(figsize=[10, 5])

plt.plot(df.dw_concs, df.index/len(df), ) 
plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
plt.xlabel('Piek drinkwaterconcentratie (g/m3)')
plt.ylabel('Cumulative kansdichtheid')
plt.title('Overschrijdingen per jaar: '+str(round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 1))+ '%, total sims:'+ str(len(df)) )
plt.xscale('log')
plt.xlim(1e-12, 1)
plt.legend()
plt.savefig(save_results_to +'/peak_monte_carlo_simulation.png', dpi=300, bbox_inches='tight')

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

df_mean = load_pickle(filename='mean_monte_carlo_loop_df', foldername='monte_carlo_output')
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

df_peak = load_pickle(filename='peak_monte_carlo_loop_df', foldername='monte_carlo_output')

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


with pd.ExcelWriter('monte_carlo_output/results_and_correlation.xlsx') as writer:  
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