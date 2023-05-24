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
                                       plot_xlog=False, 
                                       show_plot=True):
    save_results_to = check_create_folders(folder_name=folder_name)

    dfx = df.sort_values(by=parameter)
    dfx.reset_index(inplace = True, drop=True)
    fig = plt.figure(figsize=[8, 5])
    plt.plot(dfx[parameter], dfx.index/len(df), color = 'blue')
    if xlabel is None:
        plt.xlabel(parameter)
    else:
        plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if plot_xlog:
        plt.xscale('log')
    plt.savefig(save_results_to+'/' + parameter +'.png', dpi=300, bbox_inches='tight')
    if show_plot:
        plt.show()
#%%
dw_norm = 0.001

# %%
# -------------------------------------------------
# Plotting of the input parameters
# -------------------------------------------------
df = load_pickle(filename='mean_monte_carlo_F_as=3', foldername='monte_carlo_output')
save_results_to = check_create_folders(folder_name='parameter_figures')

make_plots =False

if make_plots:

    plot_distribution_input_parameter (df=df, 
                                    parameter='flow_rates', 
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative frequenteverdeling', 
                                    xlabel = '$Dagelijks huishoudelijk waterverbruik (m^{3}/dag)$',
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='activattion_energies', 
                                    xlabel= 'Activation engergie',
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='partitioning_enthalpies', 
                                    xlabel ='Partitioning enthalpie',
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Ktemps', 
                                    xlabel ="$f_{K}^{temp}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Dtemps', 
                                    xlabel ="$f_{D}^{temp}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)                                   

    plot_distribution_input_parameter (df=df, 
                                    parameter='gw_concs', 
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Kconcs', 
                                    xlabel ="$f_{K}^{conc}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='f_Dconcs', 
                                    xlabel ="$f_{D}^{conc}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)


    plot_distribution_input_parameter (df=df, 
                                    parameter='a_c_Ks', 
                                    xlabel ="$a_{KC}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='a_c_Ds', 
                                    xlabel ="$a_{DC}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)


    plot_distribution_input_parameter (df=df, 
                                    parameter='log_Kpw_refs', 
                                    xlabel ="$logK_{pw}^{ref}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='log_Dp_refs', 
                                    xlabel ="$logD_{p}^{ref}$",
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                    parameter='length_middle_points', 
                                    folder_name = save_results_to,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

    plot_distribution_input_parameter (df=df, 
                                        parameter='length_plumes', 
                                        folder_name = save_results_to,
                                        ylabel = 'Cumulative frequenteverdeling van gemeten waardes', 
                                        show_plot=True)


    # Pipe plottling
    df_PE40 = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40')

    plot_distribution_input_parameter (df=df_PE40, 
                                        parameter='Binnendiam', 
                                        xlabel='Binnen diameter (mm)',
                                        folder_name = save_results_to,
                                        ylabel = 'Cumulative kansdichtheid', 
                                        show_plot=True)

    plot_distribution_input_parameter (df=df_PE40, 
                                        parameter='Buitendiam', 
                                        xlabel='Buiten diameter (mm)',
                                        folder_name = save_results_to,
                                        ylabel = 'Cumulative kansdichtheid', 
                                        show_plot=True)

    plot_distribution_input_parameter (df=df_PE40, 
                                        parameter='Wanddikte', 
                                        xlabel='Wanddikte (mm)',
                                        folder_name = save_results_to,
                                        ylabel = 'Cumulative kansdichtheid', 
                                        show_plot=True)
#%%

plot_distribution_input_parameter (df=df_PE40, 
                                    parameter='Lengte_GIS', 
                                    xlabel='Lengte (m)',
                                    folder_name = save_results_to,
                                    plot_xlog=True,
                                    ylabel = 'Cumulative kansdichtheid', 
                                    show_plot=True)

#%%
# See that for a specific inner diameter there is a specific thickness 
df_PE40.groupby([ 'Binnendiam',]).agg(['mean', 'median', 'std', 'count', 'min', 'max' ])
df_PE40.describe()

#%% 

df_plume_concs = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')

plt.plot(df_plume_concs, np.arange(start=0, stop=len(df_plume_concs)), color = 'blue')
plt.xlabel('Plume conc (mg/kg)')
plt.ylabel('Cumulative kansdichtheid')
# plt.xscale('log')
plt.savefig(save_results_to+'/plume_conc.png', dpi=300, bbox_inches='tight')

#%%
# # range lenth plume
plume_length_values = load_pickle(filename='monte_carlo_plume_lengths', foldername='inputs')

df_plume = pd.DataFrame (plume_length_values, columns = ['plume_length_values'])
df_plume_length_values = df_plume.sort_values(by='plume_length_values')
df_plume_length_values.reset_index(inplace = True, drop=True)

plt.plot(df_plume_length_values.plume_length_values, 
         df_plume_length_values.index/len(df_plume_length_values), color = 'blue')
plt.xlabel('Plume Lengte (m)')
plt.ylabel('Cumulative verdeling van gemeten waardes')
# plt.xscale('log')
plt.savefig(save_results_to+'/lengte_plume.png', dpi=300, bbox_inches='tight')
#%%