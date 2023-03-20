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
    # size of plume (length of segment) (excel file) -> Amitosh/Auliato discuss more complicated version later
    # concentration of plume (excel file)
    # partitioning coefficient, K_ref, +/- 0.5 st. dev, @Martin to give the real values later
    # diffusion coefficient, K_ref, +/- 0.5 st. dev, @Martin to give the real values later
    # Assessment factor = 1 or 3
    # Flow rate -> @Martin, contact dw companies for distribution 
    # Pipe thickness ? -> Amitosh/Aulia

# Steps
# import the data on plume concentration, contact length
# set range for diffusion coefficient, D_ref +/- XX -> start, stop, step See table 5-1 KWR 2016.056
# set range for partitioning coefficient, K_ref +/- XX -> start, stop, step See table 5-1 KWR 2016.056
# set range flow rate
# set range pipe thickness

# Normal distrbution from mean (mu) and standard deviation (sigma), p is the 
# random number from the random number generator between 0 and 1
# NormalDist(mu, sigma).inv_cdf(p)

#%%
# import the data on plume concentration
df = pd.read_excel(module_path / 'research' / '20190702 kans normoverschrijding.xlsx', 
                    sheet_name='RGW_AH', header=[175], usecols = "A:E", ) 
ext_values = list(df.ext_value)

# range lenth pipe 
length_range =pd.read_excel(module_path / 'research' / '20190702 kans normoverschrijding.xlsx', 
                    sheet_name='RGW_AH', header=[0], usecols = "U", nrows=12) 
length_values = list(length_range.contactlengte)

#%%
# Function to create pipe
def create_pipe(ext_values, 
                length_values, 
                assessment_factor = 3):
    '''
    Create the pipe and set the conditions for the Monte-Carlo simulation
    '''

    gw_conc = random.choice(ext_values)

    # ah_todo
    # Update with information from Amitosh/Aulia
    length = random.choice(length_values)

    # set range for diffusion coefficient, for benzene
    #@Martin, is this the LogKpw_ref value update or the LogKpw?
    log_Dp_ref = -11.54717333172
    log_Dp = NormalDist(mu=log_Dp_ref, sigma=0.5).inv_cdf(p=random.random())
    
    # set range for partitioning coefficient, for benzene
    #@Martin, is this the LogDp_ref value update or the LogDp?
    log_Kpw_ref = 1.6476099999999998
    log_Kpw = NormalDist(mu=log_Kpw_ref, sigma=0.5).inv_cdf(p=random.random())
    
    #ah_todo
    # Update with information from Mirjam
    flow_rate = NormalDist(mu=0.5, sigma=0.1).inv_cdf(p=random.random())
    #ah_todo change this to be the 1,2,4 person households at 120 L per person per day

    # ah_todo
    # Update with information from Amitosh/Aulia
    wall_thickness = NormalDist(mu=0.0027, sigma=0.0001).inv_cdf(p=random.random())

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=length, #here set the length
                    inner_diameter=0.0196,
                    wall_thickness=wall_thickness)

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(concentration_groundwater = gw_conc, #here set the gw conc
                        chemical_name="Benzeen", 
                        temperature_groundwater=12,
                        flow_rate=flow_rate, 
                        suppress_print=True)

    # Set the Kpw and Dp
    seg1.log_Kpw = log_Kpw
    seg1.log_Dp= log_Dp

    #set assessment_factor
    pipe1.ASSESSMENT_FACTOR_GROUNDWATER = assessment_factor

    pipe1.validate_input_parameters()

    return pipe1, seg1

    # dw_conc = pipe1.calculate_mean_dw_concentration()
    # return gw_conc, log_Kpw, log_Dp, length, dw_conc

#%%
save_results_to = check_create_folders(folder_name='monte-carlo_output')

# Loop through the combinations, save the dw concentrations
dw_concs = []
gw_concs = []
log_Kpws = []
log_Dps = []
lengths = []
flow_rates = []
wall_thicknesses =[]

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
        pipe1, seg1 = create_pipe(ext_values=ext_values, length_values=length_values)

        dw_conc = pipe1.calculate_mean_dw_concentration()

        # gw_conc, log_Kpw, log_Dp, length, dw_conc = create_pipe(ext_values=ext_values, length_values=length_values)
        dw_concs.append(dw_conc)
        gw_concs.append(pipe1.concentration_groundwater)  
        log_Kpws.append(seg1.log_Kpw)  
        log_Dps.append(seg1.log_Dp)
        lengths.append(seg1.length)
        flow_rates.append(pipe1.flow_rate)
        wall_thicknesses.append(seg1.wall_thickness)

    # check if the 10th and 90th percentile within tolerance, then stop the loop
    tenth_perc = np.percentile(dw_concs, 10)
    ninety_perc = np.percentile(dw_concs, 90)

    criteria_ten = abs(1 - tenth_perc / tenth_perc_n_min_1)
    criteria_nine = abs(1 - ninety_perc / ninety_perc_n_min_1)

    if (criteria_ten <= tolerance) and (criteria_nine <= tolerance):
        break
    else: 
        ninety_perc_n_min_1 = ninety_perc
        tenth_perc_n_min_1 = tenth_perc

    print('ninety_perc:', ninety_perc, 'tenth_per:c', tenth_perc)

# put the data into a df, sort by the dw_conc and save
data = zip(dw_concs, gw_concs, log_Kpws, log_Dps, lengths, flow_rates, wall_thicknesses)
df = pd.DataFrame(data,  columns = ['dw_concs', 'gw_conc', 'Kpw', 'Dp', 'Length', 'flow_rates', 'wall_thicknesses'])
df.sort_values(by=['dw_concs'], inplace=True)
df.reset_index(inplace=True)

save_df_pickle(filename='example_monte-carlo_loop_df', df= df, foldername='monte-carlo_output')

#%%
save_results_to = check_create_folders(folder_name='figures')

df = load_pickle(filename='example_monte-carlo_loop_df', foldername='monte-carlo_output')
dw_norm = 0.001 # Benzene drinking water norm, g/m3

fig = plt.figure(figsize=[10, 5])

plt.plot(df.dw_concs, df.index/len(df), ) 
plt.vlines(x=dw_norm, ymin=0, ymax =1, colors='r', linestyles='--', label = 'DW Norm')
plt.xlabel('Mean DW concentration (g/m3)')
plt.ylabel('Probability')
plt.xscale('log')
plt.legend()
# plt.savefig(save_results_to +'/example_monte-carlo_simulation.png', dpi=300, bbox_inches='tight')

#%%

len(df.loc[df.dw_concs > dw_norm]) / len(df)
# %%
