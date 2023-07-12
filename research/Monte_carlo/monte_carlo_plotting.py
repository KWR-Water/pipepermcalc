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

import seaborn as sns

def plot_cumulative_distribution(df, dw_norm, save_name, save_results_to):
        fig = plt.figure(figsize=[8, 5])

        plt.plot(df.dw_concs, df.index/len(df), ) 
        plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
        plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
        plt.ylabel('Cumulative kansdichtheid')
        plt.title('Overschrijdingen per jaar: '+str(round(len(df.loc[df.dw_concs > dw_norm]) / len(df)*100, 1))+ '%, n:'+ str(len(df)) )
        plt.xscale('log')
        plt.xlim(1e-12, 1)
        plt.legend()
        plt.savefig(save_results_to +'/'+save_name + '.png', dpi=300, bbox_inches='tight')

#%%#%%
# Redo the plots



#%%
save_results_to = check_create_folders(folder_name='monte_carlo_output')

dw_norm = 0.001


def plot_cum_dist_together(simulation_name, dw_norm = 0.001):
    df_peak =load_pickle(filename='peak_monte_carlo_'+simulation_name, foldername=save_results_to)
    df_mean =load_pickle(filename='mean_monte_carlo_'+simulation_name, foldername=save_results_to)

    fig = plt.figure(figsize=[8, 5])
    plt.plot(df_mean.dw_concs, df_mean.index/len(df_mean),color = 'green', label = 'gemiddelde') 
    plt.plot(df_peak.dw_concs, df_peak.index/len(df_peak), color = 'blue', label = 'stagnatie' ) 

    plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
    plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
    plt.ylabel('Cumulative kansdichtheid')
    plt.title('Overschrijdingen per jaar: gemiddelde: '+str(round(len(df_mean.loc[df_mean.dw_concs > dw_norm]) / len(df_mean)*100, 1))+ 
            '%, stagnatie:'+str(round(len(df_peak.loc[df_peak.dw_concs > dw_norm]) / len(df_peak)*100, 1)) + '%')
    plt.xscale('log')
    plt.xlim(1e-12, 1)
    # plt.xlim([0, 0.01])
    plt.legend()
    plt.grid()
    plt.savefig(save_results_to +'/'+simulation_name + '.png', dpi=300, bbox_inches='tight')


def calculate_quantiles(simulation_name):
    df_peak =load_pickle(filename='peak_monte_carlo_'+simulation_name, foldername=save_results_to)
    df_mean =load_pickle(filename='mean_monte_carlo_'+simulation_name, foldername=save_results_to)

    peak_vals= []
    mean_vals = []

    quantile_vals = [0.01, 0.05,0.1,0.5,0.95, 0.99 ]

    for quant in quantile_vals:
        peak_vals.append(df_peak.dw_concs.quantile(quant))
        mean_vals.append(df_mean.dw_concs.quantile(quant))

        # put the data into a df, sort by the dw_conc and save
        data = zip(quantile_vals,  mean_vals, peak_vals, )
        df = pd.DataFrame(data,  
                        columns = ['quantile', 'mean', 'peak'] )
        
        df['mean/peak'] = df['mean']/df['peak']

    df.to_excel(save_results_to+ "/" +simulation_name+".xlsx")

    return df

# simulation_name = 'F_as=7'
# plot_cum_dist_together(simulation_name=simulation_name, dw_norm = 0.001)

#%%    

name_simulations = ['F_as=1', 
                     'F_as=3',
                      'F_as=7', 
                    'F_as=cum_dist_median', 
                    'F_as=3_PE80',
                     'F_as=1_median',                                          
                      ]

for simulation_name in name_simulations:

    # plot_cum_dist_together(simulation_name=simulation_name, dw_norm = 0.001)
    # calculate_quantiles(simulation_name=simulation_name,)
# %%

simulation_name = 'F_as=1'
simulation_label = 'Scenario A'


df_peak =load_pickle(filename='peak_monte_carlo_'+simulation_name, foldername=save_results_to)
df_mean =load_pickle(filename='mean_monte_carlo_'+simulation_name, foldername=save_results_to)

len_df = min(len(df_peak), len(df_mean))
df_peak = df_peak.sort_values('index')
df_peak = df_peak[1:len_df]
df_mean = df_mean.sort_values('index')
df_mean = df_mean[1:len_df]

df_mean['peak_dw'] = df_peak.dw_concs.values
df_mean['peak_gw'] = df_peak.gw_concs.values

df_mean['stagnation_factor_f'] = (((df_mean['Dp'] + 12.5) / 2 + 
                               df_mean['Kpw']) * 0.73611 + 
                                -1.03574 )

df_mean ['zero'] = 0
df_mean ["C"] = df_mean[["stagnation_factor_f", "zero"]].max(axis=1)
df_mean['stagnation_factor'] = 10 ** df_mean['C']
df_mean = df_mean.sort_values('stagnation_factor')

#%%
fig = plt.figure(figsize=[8, 5])
sns.scatterplot(data = df_mean, x='dw_concs', y = 'peak_dw',alpha = 0.5,label = simulation_label,)
# plt.scatter(df_mean.dw_concs, df_mean.peak_dw, label = simulation_label, alpha = 0.5) 
# plt.scatter(df_mean.gw_concs, df_mean.peak_gw, label = simulation_label, alpha = 0.5) 

# plt.scatter(df_mean.dw_concs, df_peak.dw_concs, label = simulation_label, alpha = 0.5, color = 'r') 
# plt.scatter(df_mean.gw_concs, df_peak.gw_concs, label = simulation_label, alpha = 0.5, color = 'r') 

plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
plt.ylabel('Stagnatie drinkwaterconcentratie (g/m3)')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim(1e-13, 1)
# plt.ylim(1e-13, 1)
plt.legend()
plt.grid()


#%%
def scatter_plot_mean_vs_peak(simulation_name, simulation_label):

    df_peak =load_pickle(filename='peak_monte_carlo_'+simulation_name, foldername=save_results_to)
    df_mean =load_pickle(filename='mean_monte_carlo_'+simulation_name, foldername=save_results_to)

    len_df = min(len(df_peak), len(df_mean))
    df_peak = df_peak.sort_values('index')
    df_peak = df_peak[1:len_df]
    df_mean = df_mean.sort_values('index')
    df_mean = df_mean[1:len_df]


    fig = plt.figure(figsize=[8, 5])
    plt.scatter(df_mean.dw_concs, df_peak.dw_concs, label = simulation_label, alpha = 0.5) 
    plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
    plt.ylabel('Stagnatie drinkwaterconcentratie (g/m3)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-13, 1)
    plt.ylim(1e-13, 1)
    plt.legend()
    plt.grid()
    plt.savefig(save_results_to +'/'+simulation_label + '.png', dpi=300, bbox_inches='tight')

scatter_plot_mean_vs_peak(simulation_name = 'F_as=1', simulation_label = 'Scenario A')
scatter_plot_mean_vs_peak(simulation_name = 'F_as=cum_dist_median', simulation_label = 'Scenario B')
scatter_plot_mean_vs_peak(simulation_name = 'F_as=3', simulation_label = 'Scenario C')
scatter_plot_mean_vs_peak(simulation_name = 'F_as=7', simulation_label = 'Scenario D')
scatter_plot_mean_vs_peak(simulation_name = 'F_as=3_PE80',simulation_label = 'Scenario E')
scatter_plot_mean_vs_peak(simulation_name = 'F_as=1_median',   simulation_label = 'Scenario F')


#%%
# def plot_cum_dist_together(simulation_name, dw_norm = 0.001):

fig = plt.figure(figsize=[8, 5])

for simulation_name in name_simulations:
    # df_peak =load_pickle(filename='peak_monte_carlo_'+simulation_name, foldername=save_results_to)
    df_mean =load_pickle(filename='mean_monte_carlo_'+simulation_name, foldername=save_results_to)

    plt.plot(df_mean.dw_concs, df_mean.index/len(df_mean), label = simulation_name) 
    # plt.plot(df_peak.dw_concs, df_peak.index/len(df_peak), color = 'blue', label = 'Peak' ) 

plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
plt.ylabel('Cumulative kansdichtheid')
plt.title('All mean')

# plt.title('Overschrijdingen per jaar: mean: '+str(round(len(df_mean.loc[df_mean.dw_concs > dw_norm]) / len(df_mean)*100, 1))+ '%')
        # '%, peak:'+str(round(len(df_peak.loc[df_peak.dw_concs > dw_norm]) / len(df_peak)*100, 1)) + '%')
plt.xscale('log')
plt.xlim(1e-12, 1)
plt.legend()
plt.grid()
plt.savefig(save_results_to +'/all_mean.png', dpi=300, bbox_inches='tight')

fig = plt.figure(figsize=[8, 5])

for simulation_name in name_simulations:
    df_peak =load_pickle(filename='peak_monte_carlo_'+simulation_name, foldername=save_results_to)
    # df_mean =load_pickle(filename='mean_monte_carlo_'+simulation_name, foldername=save_results_to)

    # plt.plot(df_mean.dw_concs, df_mean.index/len(df_mean), label = simulation_name) 
    plt.plot(df_peak.dw_concs, df_peak.index/len(df_peak), label = simulation_name) 

plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
plt.xlabel('Gemiddelde drinkwaterconcentratie (g/m3)')
plt.ylabel('Cumulative kansdichtheid')
plt.title('All peak')
# plt.title('Overschrijdingen per jaar: mean: '+str(round(len(df_mean.loc[df_mean.dw_concs > dw_norm]) / len(df_mean)*100, 1))+ '%')
        # '%, peak:'+str(round(len(df_peak.loc[df_peak.dw_concs > dw_norm]) / len(df_peak)*100, 1)) + '%')
plt.xscale('log')
plt.xlim(1e-12, 1)
plt.legend()
plt.grid()
plt.savefig(save_results_to +'/all_peak.png', dpi=300, bbox_inches='tight')