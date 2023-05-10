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
# plume_concs = list(df['Extrapolated values'])
# save_df_pickle(filename='monte_carlo_plume_concs', df= plume_concs, foldername='inputs')
plume_concs = load_pickle(filename='monte_carlo_plume_concs', foldername='inputs')

# # range lenth plume
# plume_length_range =pd.read_excel('inputs/NBNL_plume_SA.xlsx', 
#                     sheet_name='plume_SA', header=[3], usecols = "U", nrows=12) 
# plume_length_range = list(plume_length_range['Diameter (m)'])
# save_df_pickle(filename='monte_carlo_plume_lengths', df= plume_length_range, foldername='inputs')
plume_length_values = load_pickle(filename='monte_carlo_plume_lengths', foldername='inputs')

# location middle point: continuous equal distribution 0-1
length_fraction_middle_point = list(np.arange(start=0, stop=1.01, step=0.01))

# Import the pipe information - inner diameter, thickness, length
# df = pd.read_excel('Aansluitleidingen_inBedrijf_16012023 PWN.xlsx')
# save_df_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN', df=df)
# df = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN')

# df = df[['Materiaal', 'Materiaals', 
#          'Nominale_D','Buitendiam', 
#         'Binnendiam', 
#         'Wanddikte', 
#         'Lengte_GIS',  'Shape_Leng']]

# df_PE40 = df.loc[df.Materiaals == 'PE 40 (ZPE)']

# # Check that the columns are the same:
# df_PE40 .loc[df_PE40['Nominale_D'] != df_PE40['Buitendiam']]
# df_PE40['Lengte_GIS/Shape_Leng'] = df_PE40['Lengte_GIS'] / df_PE40['Shape_Leng']
# df_PE40 = df_PE40[[ 'Materiaal', 'Buitendiam', 'Binnendiam', 'Wanddikte', 'Lengte_GIS']]

# save_df_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40', df=df_PE40)

df_PE40 = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40')

# df_PE40 = load_pickle(filename='Aansluitleidingen_inBedrijf_16012023_PWN_PE40')
# df_PE40_Lengte_GIS = df_PE40.sort_values(by='Lengte_GIS')
# pipe_length_values = list(df_PE40_Lengte_GIS.Lengte_GIS)
# save_df_pickle(filename='monte_carlo_pipe_length_values', df= pipe_length_values, foldername='inputs')
pipe_length_values = load_pickle(filename='monte_carlo_pipe_length_values', foldername='inputs')

# range inner diameter
# df_PE40_inner_diam = df_PE40[df_PE40['Binnendiam'] > 0]  
# df_PE40_inner_diam = df_PE40_inner_diam.sort_values(by='Binnendiam')
# df_PE40_inner_diam['Binnendiam'] = df_PE40_inner_diam['Binnendiam']  / 1000
# df_PE40_inner_diam = df_PE40_inner_diam.round({'Binnendiam': 4})
# inner_diam_values = list(df_PE40_inner_diam.Binnendiam)
# save_df_pickle(filename='monte_carlo_df_PE40_inner_diam', df= inner_diam_values, foldername='inputs')
inner_diam_values = list(load_pickle(filename='monte_carlo_df_PE40_inner_diam', foldername='inputs'))

# dictionary relating specific inner diameters to specific wall thicknesses    
# wall_thickness_dict = {0.0124: 0.0018,
#                         0.0156: 0.0022,
#                         0.0196: 0.0027,
#                         0.0250: 0.0035,
#                         0.0314: 0.0043,
#                         0.0392: 0.0054,
#                         0.0494: 0.0068, }
# save_df_pickle(filename='monte_carlo_wall_thickness_dict', df= wall_thickness_dict, foldername='inputs')
wall_thickness_dict = load_pickle(filename='monte_carlo_wall_thickness_dict', foldername='inputs')

# Water use per household => 1/3 = 1 person, 1/3 = 2 people, 1/9 each 3,4,5 people @ 125 L/p/d
# water_use = [0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.375, 0.500, 0.625] 
# household_distribution = [1,1,1,1,1,1,1,1,1,1,
#                             2,2,2,2,2,2,2,
#                             3,3,3,
#                             4,4,4,
#                             5]

# household_water_use = 0.125 # L/p/d
# water_use =  (np.array(household_distribution) * household_water_use).tolist()
# save_df_pickle(filename='monte_carlo_water_use', df= water_use, foldername='inputs')
water_use = load_pickle(filename='monte_carlo_water_use', foldername='inputs')


input_parameters = {
                    'concentration_soil': plume_concs, 
                    'length_pipe':  pipe_length_values,
                    'length_fraction_middle_point': length_fraction_middle_point,
                    'length_plume':  plume_length_values,
                    'inner_diameter': inner_diam_values, 
                    'flow_rate': water_use, 
                    'wall_thickness_dict': wall_thickness_dict,
                    }

save_df_pickle(filename='monte_carlo_input_parameters', df= input_parameters, foldername='inputs')

#%%