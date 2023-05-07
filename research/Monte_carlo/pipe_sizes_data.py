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

#%%
# -------------------------------------------------
# Plotting
# -------------------------------------------------

save_results_to = check_create_folders(folder_name='figures')

# Binnen diameter (mm)
dfx = df_PE40.sort_values(by='Binnendiam')
dfx.reset_index(inplace = True, drop=True)
plt.plot(dfx.Binnendiam, dfx.index/len(dfx), color = 'blue')
plt.xlabel('Binnen diameter (mm)')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/Binnen_diameter.png', dpi=300, bbox_inches='tight')

#%%
# Buiten diameter (mm)
dfx = df_PE40.sort_values(by='Buitendiam')
dfx.reset_index(inplace = True, drop=True)
plt.plot(dfx.Buitendiam, dfx.index/len(dfx), color = 'blue')
plt.xlabel('Buiten diameter (mm)')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/Buiten_diameter.png', dpi=300, bbox_inches='tight')

#%%
# Wanddikte (mm)
df_PE40_Wanddikte = df_PE40.sort_values(by='Wanddikte')
df_PE40_Wanddikte.reset_index(inplace = True, drop=True)
plt.plot(df_PE40_Wanddikte.Wanddikte,dfx.index/len(dfx), color = 'blue')
plt.xlabel('Wanddikte (mm)')
plt.ylabel('Cumulatieve dichtheid')
plt.savefig(save_results_to+'/wanddikte.png', dpi=300, bbox_inches='tight')

#%%
# Lengte_GIS (mm)
df_PE40_Lengte_GIS = df_PE40.sort_values(by='Lengte_GIS')
df_PE40_Lengte_GIS.reset_index(inplace = True, drop=True)
plt.plot(df_PE40_Lengte_GIS.Lengte_GIS, dfx.index/len(dfx), color = 'blue')
plt.xlabel('Lengte (m)')
plt.ylabel('Cumulatieve dichtheid')
plt.xscale('log')
plt.savefig(save_results_to+'/lengte.png', dpi=300, bbox_inches='tight')

#%%
# See that for a specific inner diameter there is a specific thickness 
df_PE40.groupby([ 'Binnendiam',]).agg(['mean', 'median', 'std', 'count', 'min', 'max' ])
df_PE40.describe()

#%%



#%%