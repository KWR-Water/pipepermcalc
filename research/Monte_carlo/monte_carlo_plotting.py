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

def plot_distribution_input_parameter (df, 
                                       parameter, 
                                       folder_name, 
                                       xlabel=None,
                                       ylabel = 'Cumulative kansdichtheid', 
                                       show_plot=True):
    save_results_to = check_create_folders(folder_name=folder_name)

    dfx = df.sort_values(by=parameter)
    dfx.reset_index(inplace = True, drop=True)
    fig = plt.figure(figsize=[10, 5])
    plt.plot(dfx[parameter], dfx.index/len(df), color = 'blue')
    if xlabel is None:
        plt.xlabel(parameter)
    else:
        plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_results_to+'/' + parameter +'.png', dpi=300, bbox_inches='tight')
    if show_plot:
        plt.show()
#%%
dw_norm = 0.001

# %%
# -------------------------------------------------
# Plotting of the input parameters
# -------------------------------------------------
df = load_pickle(filename='mean_monte_carlo_loop_df', foldername='monte_carlo_output')
save_results_to = check_create_folders(folder_name='figures')

make_plots =False

if make_plots:

    plot_distribution_input_parameter (df=df, 
                                    parameter='flow_rates', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative frequenteverdeling', 
                                    xlabel = 'Dagelijks huishoudelijk waterverbruik (m3/dag)',
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='activattion_energies', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='partitioning_enthalpies', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Ktemps', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Dtemps', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)                                   

    plot_distribution_input_parameter (df=df, 
                                    parameter='gw_concs', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Kconcs', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Dconcs', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)


    plot_distribution_input_parameter (df=df, 
                                    parameter='a_c_Ks', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='a_c_Ds', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)


    plot_distribution_input_parameter (df=df, 
                                    parameter='log_Kpw_refs', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='log_Dp_refs', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='length_middle_points', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)


plot_distribution_input_parameter (df=df, 
                                    parameter='length_plumes', 
                                    folder_name='figures', 
                                    ylabel = 'Cumulative frequenteverdeling van gemeten waardes', 
                                    show_plot=True)

#%% ###########################

df_plume_concs = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')
df_plume_concs = df_plume_concs.sort_values(by='Extrapolated values')
df_plume_concs.reset_index(inplace = True, drop=True)

plt.plot(df_plume_concs['Extrapolated values'], df_plume_concs.index/len(df_plume_concs), color = 'blue')
plt.xlabel('Plume conc (mg/kg)')
plt.ylabel('Cumulative kansdichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/plume_conc.png', dpi=300, bbox_inches='tight')

#%%
# range lenth plume
# plume_length_range =pd.read_excel('inputs/20190702 kans normoverschrijding.xlsx', 
#                     sheet_name='RGW_AH', header=[0], usecols = "U", nrows=12) 
# save_df_pickle(filename='monte_carlo_plume_lengths', df= plume_length_range, foldername='inputs')
plume_length_range = load_pickle(filename='monte_carlo_plume_lengths', foldername='inputs')
plume_length_values = list(plume_length_range['Diameter (m)'])

df_plume = pd.DataFrame (plume_length_values, columns = ['plume_length_values'])
df_plume_length_values = df_plume.sort_values(by='plume_length_values')
df_plume_length_values.reset_index(inplace = True, drop=True)

plt.plot(df_plume_length_values.plume_length_values, 
         df_plume_length_values.index/len(df_plume_length_values), color = 'blue')
plt.xlabel('Plume Lengte (m)')
plt.ylabel('Cumulative verdeling van gemeten waardes')
# plt.xscale('log')
plt.savefig(save_results_to+'/lengte_plume.png', dpi=300, bbox_inches='tight')