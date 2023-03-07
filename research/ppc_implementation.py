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

#%%
seg1 = Segment(name='seg1',
            material= 'PE40',
            length=25,
            inner_diameter=0.0196,
            wall_thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1])
input_gw = 1

df = read_csv('did_not_pass_chems_peak.csv')
failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

database = pipe1.view_database_chemical_names( language='NL')
database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])
database = database.loc[database['log_distribution_coefficient']>=0]
database_chemicals = database['chemical_name_NL']

# database_chemicals= ['benzeen']
fails = []
for chemical_name in database_chemicals:

# chemical_name= 'benzo[b]fluoranthene'



    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    counter = 0
    max_iterations = 1000
    relaxation_factor = 0.5
    tolerance = 0.01
    concentration_groundwater = input_gw
    concentration_drinking_water_n_min_1 = 0
    concentration_drinking_water_n_plus_1 = 0
    lower_limit = 0
    input_gw = 1
    upper_limit = input_gw
    criteria_list = [0]
    min_criteria = 100
    proc = 0

    while True:    
        concentration_drinking_water_n_min_1 = concentration_drinking_water_n_plus_1
        concentration_groundwater = input_gw

        sum_mass_segment = 0

        for segment in pipe1.segment_list:
            segment._calculate_peak_dw_mass_per_segment(pipe=pipe1, 
                                                        concentration_drinking_water=concentration_drinking_water_n_min_1,
                                    concentration_groundwater=pipe1.concentration_groundwater,)

            sum_mass_segment += segment.mass_chemical_drinkwater

        concentration_drinking_water_n = (sum_mass_segment / 
                                        pipe1.total_volume ) 
        counter +=1

        criteria = abs(1 - concentration_drinking_water_n_min_1 / concentration_drinking_water_n) / relaxation_factor

        if criteria <= tolerance:
            print('solution_found!')
            break
        elif counter > max_iterations:
            fails.append(chemical_name)
            print('Max iterations exceeded')
            break
        else:
            min_criteria = min(min_criteria, criteria)
            criteria_list.append(criteria)

            if counter == 1:
                concentration_drinking_water_n_plus_1 = concentration_groundwater *0.999
            if counter == 2:
                concentration_drinking_water_n_plus_1 = concentration_groundwater * 0.0001
            if counter >2:
                if (criteria < criteria_list[counter-1]) or (concentration_drinking_water_n > concentration_groundwater):
                    lower_limit = concentration_drinking_water_n_min_1
                    concentration_drinking_water_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                else:
                    upper_limit = concentration_drinking_water_n_min_1
                    concentration_drinking_water_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                    
            # print(concentration_drinking_water_n_min_1, concentration_drinking_water_n, criteria, proc, lower_limit, upper_limit)
            criteria_old = criteria

        # if counter % 100 ==0 : print(concentration_drinking_water_n_min_1, concentration_drinking_water_n)

len(fails)


#%%
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
database_chemicals = database['chemical_name_NL']

df = read_csv('did_not_pass_chems_peak.csv')
failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

for chemical_name in failed_chem:
    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    mean_conc=pipe1.calculate_peak_dw_concentration(relaxation_factor = 0.7, 
                                                    tolerance = 0.001)


    pipe1.set_conditions(chemical_name=chemical_name, 
                        temperature_groundwater=12, 
                        concentration_drinking_water = mean_conc,
                        flow_rate=0.5)

    output_gw = pipe1.calculate_peak_allowable_gw_concentration(tolerance = 0.001)

    if abs(1-(input_gw/output_gw)) < 0.001:
        chem_name.append(chemical_name)
        pass_fail.append('passed')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append(output_gw)
    else: 
        chem_name.append(chemical_name)
        pass_fail.append('failed')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append(output_gw)

df = pd.DataFrame(list(zip(chem_name, pass_fail, mean_conc_vals, output_gw_vals)),
               columns =['chemical_name_NL', 'pass_fail', 'mean_conc', 'output_gw'])

database_dict = database.set_index('chemical_name_NL').to_dict()

for key, value in database_dict.items():
    df[key] = df['chemical_name_NL'].map(database_dict[key])

df.to_csv('did_not_pass_chems_peak.csv')

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
pipe1.calculate_peak_dw_concentration(relaxation_factor=0.7)
#%%
#%%

concentration_drinking_water_n_plus_1 = 0
lower_limit = 0
input_gw = 1
upper_limit = input_gw
criteria_old = 0
counter = 0
#%%
chemical_name = 'antraceen'

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_conditions(
    chemical_name=chemical_name, 
                    concentration_groundwater =input_gw,
                    temperature_groundwater=12, 
                    flow_rate=0.5)

pipe1.validate_input_parameters()
# counter = 0
max_iterations = 100
relaxation_factor = 0.5
tolerance = 0.01
# concentration_drinking_water = 0.984375 #0.977
concentration_drinking_water_n_min_1 = 0.006
concentration_groundwater = input_gw

sum_mass_segment = 0

for segment in pipe1.segment_list:
    segment.stagnation_factor = segment._calculate_stagnation_factor()
    # delta_c = concentration_groundwater - concentration_drinking_water
    delta_c = concentration_groundwater - concentration_drinking_water_n_min_1

    # From equation 4-10 KWR 2016.056, but not simplifying the mass flux 
    # in equation 4-5 and rearranging to remove C_dw from the equation       
    segment.mass_chemical_drinkwater = (((10 ** segment.log_Dp * 10 ** segment.log_Kpw)
                                            * segment.permeation_surface_area 
                                            * delta_c / segment.diffusion_path_length 
                                            * pipe1.stagnation_time * segment.stagnation_factor) 
                                        / segment.ASSESSMENT_FACTOR_GROUNDWATER)

    sum_mass_segment += segment.mass_chemical_drinkwater

# concentration_pipe_drinking_water = (sum_mass_segment / 
#                                 pipe1.total_volume ) 
concentration_drinking_water_n = (sum_mass_segment / 
                                pipe1.total_volume ) 

# criteria = abs(1 - concentration_drinking_water / concentration_pipe_drinking_water)
criteria = abs(1 - concentration_drinking_water_n_min_1 / concentration_drinking_water_n)

# concentration_drinking_water, concentration_pipe_drinking_water, criteria
counter +=1
concentration_drinking_water_n_min_1, concentration_drinking_water_n, lower_limit, criteria, criteria_old, counter

#%%
if counter == 1:
    concentration_drinking_water_n_plus_1 = concentration_groundwater / 2
# elif concentration_drinking_water_n > concentration_groundwater:
elif (criteria < criteria_old) or (concentration_drinking_water_n > concentration_groundwater):
    concentration_drinking_water_n_plus_1 = concentration_drinking_water_n_min_1 + (upper_limit- lower_limit) / 2
    lower_limit = concentration_drinking_water_n_min_1
else:
    concentration_drinking_water_n_plus_1 = concentration_drinking_water_n_min_1 - (upper_limit - lower_limit ) / 2
    upper_limit = concentration_drinking_water_n_min_1
criteria_old = criteria

if criteria < tolerance:
    print('done!')

concentration_drinking_water_n_plus_1

#%%
mean_conc=pipe1.calculate_peak_dw_concentration(tolerance = 0.01, relaxation_factor = 0.9, max_iterations=1000)

if np.isnan(mean_conc):
    pass
else:
    pipe1.set_conditions(chemical_name=chemical_name, 
                        temperature_groundwater=12, 
                        concentration_drinking_water = mean_conc,
                        flow_rate=0.5)

    output_gw = pipe1.calculate_peak_allowable_gw_concentration(tolerance = 0.001)

    if abs(1-(input_gw/output_gw)) < 0.001:
        print(output_gw)
    else: 
        print(output_gw)
#%%------------------------------------------------------------------------------------    
# test mean
seg1 = Segment(name='seg1',
            material= 'PE40',
            length=25,
            inner_diameter=0.0196,
            wall_thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1])

input_gw = 10

pipe1.set_conditions(
    chemical_name="Benzeen", 
                    concentration_groundwater =input_gw,
                    temperature_groundwater=12, 
                    flow_rate=0.5 

                    )

pipe1.validate_input_parameters()

mean_conc=pipe1.calculate_mean_dw_concentration(tolerance = 0.00001)


pipe1.set_conditions(chemical_name="Benzeen", 
                    temperature_groundwater=12, 
                    concentration_drinking_water = mean_conc,
                    flow_rate=0.5)

output_gw = pipe1.calculate_mean_allowable_gw_concentration(tolerance = 0.00001)

if abs(1-(input_gw/output_gw)) < 0.001:
    print('passed')
else: print(input_gw, output_gw, input_gw/output_gw)
