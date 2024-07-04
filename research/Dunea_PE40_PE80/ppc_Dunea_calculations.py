#%% ----------------------------------------------------------------------------
# A. Hockin, January 2023
# KWR 403230-003
# Pipe permeation calculator
# With Martin vd Schans, Bram Hillebrand, Lennart Brokx
#
# ------------------------------------------------------------------------------

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from pandas import read_csv
from pandas import read_excel
from datetime import timedelta
from scipy.optimize import minimize
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

from project_path import file_path

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

# check_create_folders(folder_name='output')

def check_soil_gw_concentrations (save_name, material):

    # Make segment from the median values in KWR 2024.060
    seg1 = Segment(name='seg1',
                material= material,
                length=6.91,
                inner_diameter=0.0196,
                wall_thickness=0.0027, 
                )

    pipe1 = Pipe(segment_list=[seg1])

    # load database, keep only chemicals with a knwon drinking water norm and 
    # where the drinking water norm is less than the solubility
    database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])

    mask1 = (database['Drinking_water_norm'] < database['solubility'])
    database = database[mask1]
    database_chemicals = database['chemical_name_NL']
    drinking_norms = database['Drinking_water_norm']
    solubilities = database['solubility']

    k_coefficients = []
    d_coefficients = []
    soil_concentrations = []
    gw_concentrations = []

    log_Dp_refs =[]
    log_Kpw_refs =[]

    a_c_Ds = []
    a_c_Ks = []
    f_Dconcs =[]
    f_Kconcs =[]

    activattion_energies = []
    partitioning_enthalpies = []
    f_Dtemps =[]
    f_Ktemps =[]


    for chemical_name, drinking_norm, solubility in zip(database_chemicals, drinking_norms, solubilities):

        pipe1.set_conditions(chemical_name=chemical_name, 
                            temperature_groundwater=12, 
                            concentration_drinking_water = drinking_norm,
                            flow_rate=0.25, #median household water use, m3/day
                            suppress_print = True
                            ) 
        
        pipe1.validate_input_parameters()
        gw_concentration = pipe1.calculate_mean_allowable_gw_concentration()
        if round(gw_concentration, 8) >= solubility:
            gw_concentrations.append(np.nan)
            k_coefficients.append(np.nan)
            d_coefficients.append(np.nan)
            soil_concentrations.append(np.nan)
            log_Kpw_refs.append(np.nan)
            a_c_Ks.append(np.nan)
            f_Kconcs.append(np.nan)
            partitioning_enthalpies.append(np.nan)
            f_Ktemps.append(np.nan)           
        else:
            gw_concentrations.append(gw_concentration)
            k_coefficients.append(seg1.log_Kpw)
            d_coefficients.append(seg1.log_Dp)
            soil_concentrations.append(pipe1.concentration_soil)

            log_Kpw_refs.append(seg1.log_Kpw_ref)
            a_c_Ks.append(seg1.PARTITIONING_A_C)
            f_Kconcs.append(seg1.f_Kconc)
            partitioning_enthalpies.append(seg1.partitioning_enthalpie)
            f_Ktemps.append(seg1.f_Ktemp)

            # log_Dp_refs.append(log_Dp_ref)
            # a_c_Ds.append(a_c_D)
            # f_Dconcs.append(f_Dconc)
            # activattion_energies.append(seg1.activattion_energy)
            # f_Dtemps.append(f_Dtemp)            


    df_save= pd.DataFrame(
        {'chemical_name_NL': database_chemicals,
        'Drink water norm (g/m3)': drinking_norms,
        'Log Partitiecoëfficiënt (-)': k_coefficients,
        'Log Diffusiecoëfficiënt (m2/s)': d_coefficients,
        'Bodem gehalte (mg/kg)': soil_concentrations, 
        'Grondwater concentratie (g/m3)': gw_concentrations,
        # 'log_Dp_refs', 
        'log_Kpw_refs': log_Kpw_refs,
        # 'f_Dconcs', 
        'f_Kconcs':f_Kconcs, 
        # 'a_c_Ds', 
        'a_c_Ks':a_c_Ks,
        # 'f_Dtemps', 
        'f_Ktemps':f_Ktemps, 
        # 'activattion_energies', 
        'partitioning_enthalpies':partitioning_enthalpies
        })

    df_save.to_excel(save_name + '.xlsx')
    save_df_pickle(filename=save_name, df=df_save, foldername='')
    return df_save, database

df_pe40, database = check_soil_gw_concentrations (save_name = 'pe40', material= 'PE40')
df_pe80, database = check_soil_gw_concentrations (save_name = 'pe80', material= 'PE80')

#%%
df_save = pd.merge(df_pe40.dropna(), df_pe80.dropna(), on=["chemical_name_NL",'Drink water norm (g/m3)'] , suffixes= ('_PE40', '_PE80'))
df_save = pd.merge(database, df_save, on=["chemical_name_NL"] )
df_save['Cgw PE80/PE40'] = df_save['Grondwater concentratie (g/m3)_PE80'] / df_save['Grondwater concentratie (g/m3)_PE40']
df_save['logKpw PE80/PE40'] = df_save['Log Partitiecoëfficiënt (-)_PE80'] / df_save['Log Partitiecoëfficiënt (-)_PE40']
df_save['Log Dp PE80/PE40'] = df_save['Log Diffusiecoëfficiënt (m2/s)_PE80'] / df_save['Log Diffusiecoëfficiënt (m2/s)_PE40']

df_save.to_excel('pe40_vs_pe80.xlsx')
save_df_pickle(filename='pe40_vs_pe80', df =df_save, foldername='')

df_export = df_save[['CAS_number', 'chemical_name_NL', 'chemical_name_EN',
       'molecular_weight', 'solubility',
       'log_octanol_water_partitioning_coefficient',
       'log_distribution_coefficient', 'chemical_group',
       'chemical_group_number', 'molecular_volume', 
       'Drink water norm (g/m3)', 'Log Partitiecoëfficiënt (-)_PE40',
       'Log Diffusiecoëfficiënt (m2/s)_PE40', 'Bodem gehalte (mg/kg)_PE40',
       'Grondwater concentratie (g/m3)_PE40',  'Log Partitiecoëfficiënt (-)_PE80',
       'Log Diffusiecoëfficiënt (m2/s)_PE80', 'Bodem gehalte (mg/kg)_PE80',
       'Grondwater concentratie (g/m3)_PE80', 'Cgw PE80/PE40', 'logKpw PE80/PE40',
       'Log Dp PE80/PE40', ]]
df_export.to_excel('comparison_pe40_vs_pe80.xlsx')
save_df_pickle(filename='comparison_pe40_vs_pe80', df =df_save, foldername='')

#%%
df_save.groupby(['chemical_group_number'])['Cgw PE80/PE40','logKpw PE80/PE40', 'Log Dp PE80/PE40' ].mean()

df_save['Cgw PE80/PE40'].mean(), df_save['Cgw PE80/PE40'].std()
#%%
# For plotting
mask1 = (df_save['Cgw PE80/PE40'] <= 1)
df_save[mask1]

df_save['style'] = 1
df_save['style'][mask1] = 2
#%%
param = 'Grondwater concentratie (g/m3)'
xparameter = param + '_PE40'
yparameter = 'Cgw PE80/PE40' #param + '_PE80'
hue_param = 'chemical_group_number' #'log_octanol_water_partitioning_coefficient'
style_param = 'style'
line_slope = 1
save_name = 'PE40_v_PE80'

fig = plt.figure(figsize=[6, 6])
# plt.gca().set_aspect('equal')

sns.scatterplot(data = df_save, x =xparameter, y =yparameter, 
                style = style_param, 
                hue= hue_param,
                # palette='viridis',
                s=50,
                zorder=10)
# add a line for the relationship between the x and y parameters (if known)
plt.axline((0,0), slope=line_slope, color='black', linestyle='dashed', label='1:1')
# plt.xlim(ylims)
# plt.ylim(ylims)
plt.xlabel(xparameter)
plt.ylabel(yparameter)
# plt.yscale('log')
# plt.xscale('log')

# plt.title(df_compare_dict['sample_description'][evides_code])
plt.grid()
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)

plt.savefig('crossplot_' + str(save_name) +'.png', format='png',
                            bbox_inches='tight',
                            dpi = (500))# %%

#%%
param = 'Grondwater concentratie (g/m3)'
xparameter = param + '_PE40'
yparameter = 'Cgw PE80/PE40' #param + '_PE80'
hue_param = 'chemical_group_number' #'log_octanol_water_partitioning_coefficient'
style_param = 'style'
line_slope = 1
save_name = 'PE40_v_PE80_factor'

fig = plt.figure(figsize=[6, 6])
# plt.gca().set_aspect('equal')

sns.scatterplot(data = df_save, x =xparameter, y =yparameter, 
                # style = style_param, 
                # hue= hue_param,
                # palette='viridis',
                s=50,
                zorder=10)

mask1 = (df_save['chemical_name_NL'] == 'Benzeen')
df_plot = df_save[mask1]
plt.scatter(df_plot[xparameter], df_plot[yparameter], label = 'Benzeen', c= 'r', zorder = 20)

plt.axhline(y= 1, color='black', linestyle='dashed', label='1')
# plt.xlim(ylims)
# plt.ylim(ylims)
plt.xlabel(xparameter)
plt.ylabel('Factor (Grondwater Concentratie PE80/PE40)') #yparameter)
# plt.yscale('log')
plt.xscale('log')

# plt.title(df_compare_dict['sample_description'][evides_code])
plt.grid()
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)

plt.savefig('crossplot_' + str(save_name) +'.png', format='png',
                            bbox_inches='tight',
                            dpi = (500))# %%
#%%
param = 'Log Diffusiecoëfficiënt (m2/s)'
xparameter = param + '_PE40' #'Diffusiecoëfficiënt (m2/s)_PE40'# 
yparameter = param + '_PE80' #'Partitiecoëfficiënt (-)_PE40' 
hue_param =  'chemical_group_number' #'log_octanol_water_partitioning_coefficient' #
style_param = 'style'
line_slope = 1
save_name = 'PE40_v_PE80_diffusion'

fig = plt.figure(figsize=[6, 6])
# plt.gca().set_aspect('equal')

sns.scatterplot(data = df_save, x =xparameter, y =yparameter, 
                style = style_param, 
                hue= hue_param,
                # hue_norm=(0.002, 1),
                # palette='viridis',
                s=50,
                zorder=10)
# add a line for the relationship between the x and y parameters (if known)
plt.axline((0,0), slope=line_slope, color='black', linestyle='dashed', label='1:1')
plt.xlim([-16, -12])
plt.ylim([-16, -12])
plt.xlabel(xparameter)
plt.ylabel(yparameter)
# plt.yscale('log')
# plt.xscale('log')

# plt.title(df_compare_dict['sample_description'][evides_code])
plt.grid()
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)

plt.savefig('crossplot_' + str(save_name) +'.png', format='png',
                            bbox_inches='tight',
                            dpi = (500))# %%

#%%
# Make segment from the median values in KWR 2024.060
seg1 = Segment(name='seg1',
            material= 'PE40',
            length=6.91,
            inner_diameter=0.0196,
            wall_thickness=0.0027, 
            )

pipe1 = Pipe(segment_list=[seg1])

# load database, keep only chemicals with a known drinking water norm and 
# where the drinking water norm is less than the solubility
database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])

mask1 = (database['Drinking_water_norm'] < database['solubility'])
database = database[mask1]

pipe1.set_conditions(chemical_name='benzo(g,h,i)perylene', 
                    temperature_groundwater=12, 
                    concentration_drinking_water = 0.0001,
                    flow_rate=0.25, #median household water use, m3/day
                    # suppress_print = True
                    ) 

pipe1.validate_input_parameters()
gw_concentration = pipe1.calculate_mean_allowable_gw_concentration()

# seg1 = Segment(name='seg1',
#             material= 'PE80',
#             length=6.91,
#             inner_diameter=0.0196,
#             wall_thickness=0.0027, 
#             )

# pipe1 = Pipe(segment_list=[seg1])
# pipe1.set_conditions(chemical_name='benzo(g,h,i)perylene', 
#                     temperature_groundwater=12, 
#                     concentration_drinking_water = 0.0001,
#                     flow_rate=0.25, #median household water use, m3/day
#                     # suppress_print = True
#                     ) 

# pipe1.validate_input_parameters()
# gw_concentration = pipe1.calculate_mean_allowable_gw_concentration()

# # solubility= database['solubility'].loc[database['chemical_name_NL'] == 'benzo(g,h,i)perylene']

# # round(gw_concentration, 6) >= solubility
#%%