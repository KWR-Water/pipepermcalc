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

from project_path import file_path

from pipepermcalc.pipe import * 
from pipepermcalc.segment import * 


# chem_name = []
# pass_fail = []
# mean_conc_vals = []
# output_gw_vals = []
# passed = []
# seg1 = Segment(name='seg1',
#             material= 'PE40',
#             length=25,
#             inner_diameter=0.0196,
#             wall_thickness=0.0027,
#             )

# pipe1 = Pipe(segment_list=[seg1])
    
# pipe1.set_conditions(chemical_name='chloordaan', 
#                 temperature_groundwater=12, 
#                 concentration_drinking_water = 0.0072075941,
#                 flow_rate=0.5)

# pipe1.validate_input_parameters()

# pipe1.calculate_mean_allowable_gw_concentration(tolerance = 0.01, debug = True)

#%% PEAK
chem_name = []
pass_fail = []
mean_conc_vals = []
output_gw_vals = []
passed = []
seg1 = Segment(name='seg1',
            material= 'PE40',
            length=25,
            inner_diameter=0.0196,
            wall_thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1])
input_gw = 1

database = pipe1.view_database_chemical_names( language='NL')
database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])
database = database.loc[database['log_distribution_coefficient']>=0]
database = database.loc[database['Drinking_water_norm'] < database['solubility'] ]
database_chemicals = database['chemical_name_NL']
solubilities = database['solubility']

len(database_chemicals)
over_solubility = []
input_gw_list = []
failed = []

for chemical_name, solubiliy in zip(database_chemicals, solubilities):
    if input_gw > solubiliy:
        input_gw = 0.01 * solubiliy

    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    mean_conc=pipe1.calculate_peak_dw_concentration()

    if mean_conc > pipe1.solubility:
        input_gw_list.append(input_gw)
        chem_name.append(chemical_name)
        pass_fail.append('over_solubilty')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append('xx')
    else:
        
        pipe1.set_conditions(chemical_name=chemical_name, 
                            temperature_groundwater=12, 
                            concentration_drinking_water = mean_conc,
                            flow_rate=0.5)

        output_gw = pipe1.calculate_peak_allowable_gw_concentration(max_iterations=500)

        if abs(1-(input_gw/output_gw)) < 0.01:
            input_gw_list.append(input_gw)
            chem_name.append(chemical_name)
            pass_fail.append('passed')
            mean_conc_vals.append(mean_conc)
            output_gw_vals.append(output_gw)
        else: 
            failed.append(chemical_name)
            input_gw_list.append(input_gw)
            chem_name.append(chemical_name)
            pass_fail.append('failed')
            mean_conc_vals.append(mean_conc)
            output_gw_vals.append(output_gw)

df_peak = pd.DataFrame(list(zip(chem_name, pass_fail, input_gw_list, mean_conc_vals, output_gw_vals)),
               columns =['chemical_name_NL', 'pass_fail', 'input_gw', 'mean_conc', 'output_gw'])

database_dict = database.set_index('chemical_name_NL').to_dict()

for key, value in database_dict.items():
    df_peak[key] = df_peak['chemical_name_NL'].map(database_dict[key])

df_peak.to_csv('did_not_pass_chems_peak_allowable.csv')

#%%
# MEAN
chem_name = []
pass_fail = []
mean_conc_vals = []
output_gw_vals = []
passed = []
seg1 = Segment(name='seg1',
            material= 'PE40',
            length=25,
            inner_diameter=0.0196,
            wall_thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1])
input_gw = 1

database = pipe1.view_database_chemical_names( language='NL')
database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])
database = database.loc[database['log_distribution_coefficient']>=0]
database = database.loc[database['Drinking_water_norm'] < database['solubility'] ]
database_chemicals = database['chemical_name_NL']
solubilities = database['solubility']

failed = []
# df = read_csv('did_not_pass_chems_mean_allowable.csv')
# failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

for chemical_name, solubiliy in zip(database_chemicals, solubilities):
    if input_gw > solubiliy:
        input_gw = 0.01 * solubiliy

    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    mean_conc=pipe1.calculate_mean_dw_concentration()

    pipe1.set_conditions(chemical_name=chemical_name, 
                        temperature_groundwater=12, 
                        concentration_drinking_water = mean_conc,
                        flow_rate=0.5)
    if mean_conc > pipe1.solubility:
        chem_name.append(chemical_name)
        pass_fail.append('over_solubilty')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append('xx')
    else:
        
        pipe1.set_conditions(chemical_name=chemical_name, 
                            temperature_groundwater=12, 
                            concentration_drinking_water = mean_conc,
                            flow_rate=0.5)

        output_gw = pipe1.calculate_mean_allowable_gw_concentration(scale_factor_upper_limit= 0.9)

        if abs(1-(input_gw/output_gw)) < 0.01:
            chem_name.append(chemical_name)
            pass_fail.append('passed')
            mean_conc_vals.append(mean_conc)
            output_gw_vals.append(output_gw)
        else: 
            failed.append(chemical_name)
            chem_name.append(chemical_name)
            pass_fail.append('failed')
            mean_conc_vals.append(mean_conc)
            output_gw_vals.append(output_gw)

df_mean = pd.DataFrame(list(zip(chem_name, pass_fail, mean_conc_vals, output_gw_vals)),
               columns =['chemical_name_NL', 'pass_fail', 'mean_conc', 'output_gw'])

database_dict = database.set_index('chemical_name_NL').to_dict()

for key, value in database_dict.items():
    df_mean[key] = df_mean['chemical_name_NL'].map(database_dict[key])

df_mean.to_csv('did_not_pass_chems_mean_allowable.csv')

#%%
chemical_name = 'Fenanthreen'
chemical_name = 'antraceen'

seg1 = Segment(name='seg1',
            material= 'PE40',
            length=25,
            inner_diameter=0.0196,
            wall_thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1])
input_gw = 1

pipe1.set_conditions(
    chemical_name=chemical_name, 
                    concentration_groundwater =input_gw,
                    temperature_groundwater=12, 
                    flow_rate=0.5)

pipe1.validate_input_parameters()
pipe1.calculate_peak_dw_concentration()
#%%
#%%

