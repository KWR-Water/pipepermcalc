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
# working version DW->GW MEAN
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
# database = database.loc[database['log_distribution_coefficient']>=0]
database_chemicals = database['chemical_name_NL']

df = read_csv('did_not_pass_chems_mean_allowable.csv')
failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

# database_chemicals= ['Hexachloorbenzeen']
# chemical_name= 'benzo[b]fluoranthene'

fails = []
over_solubility=[]
for chemical_name in database_chemicals: #database_chemicals: # new_fails: #failed_chem: # new_fails: #

    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        # concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    pipe1._fetch_chemical_database(chemical_name=pipe1.chemical_name, 
                                          suppress_print = True, 
                                          language=pipe1.language)

    if pipe1.concentration_drinking_water > pipe1.solubility:
        # raise ValueError(
            over_solubility.append(chemical_name)

            print('Error, the drinking water concentration given or the default drinking water norm is higher than the solubility of the chemical. \
            Input a lower drinking water concentration using .set_conditions()')
    else:
        # calculate initial guess for gw concentration
        sum_KDA_d = 0
        for segment in pipe1.segment_list:
            # calculate the sum of the Kpw * DP * SA *f_stag / d for all pipe segments
            log_Dp_ref = segment._calculate_ref_logD(chemical_group_number=pipe1.chemical_group_number,
                        molecular_weight=pipe1.molecular_weight)
            log_Kpw_ref = segment._calculate_ref_logK(chemical_group_number=pipe1.chemical_group_number,
                        log_octanol_water_partitioning_coefficient=pipe1.log_octanol_water_partitioning_coefficient)
            
            sum_KDA_d_segment = ( 10 ** log_Dp_ref * 10 ** log_Kpw_ref * segment.permeation_surface_area 
                                    / segment.diffusion_path_length )

            sum_KDA_d += sum_KDA_d_segment

        # initial guess concentration in groundwater
        concentration_groundwater_n_plus_1 = (pipe1.concentration_drinking_water * (1
                                         + pipe1.flow_rate * segment.ASSESSMENT_FACTOR_GROUNDWATER ) 
                                            / sum_KDA_d ) * 24* 60 * 60
        counter = 0
        max_iterations = 100
        relaxation_factor = 0.5
        tolerance = 0.01
        concentration_drinking_water_n_min_1 = pipe1.concentration_drinking_water
        concentration_drinking_water_n_plus_1 = pipe1.concentration_drinking_water
        lower_limit = pipe1.concentration_drinking_water
        upper_limit = pipe1.solubility
        criteria_list = [0]
        min_criteria = 100
        proc = 0

        while True:  
            concentration_groundwater_n_min_1 = concentration_groundwater_n_plus_1

            pipe1.set_conditions(chemical_name=pipe1.chemical_name,                                    
                concentration_groundwater=concentration_groundwater_n_min_1,
                flow_rate=pipe1.flow_rate,
                concentration_drinking_water=pipe1.concentration_drinking_water, 
                temperature_groundwater=pipe1.temperature_groundwater, 
                stagnation_time = pipe1.stagnation_time,
                suppress_print = True, 
                language = pipe1.language)

            sum_mass_segment = 0

            # mass of chemical in pipe water to meet drinking water norm
            mass_drinkingwater_norm = (pipe1.concentration_drinking_water * pipe1.flow_rate)

            for segment in pipe1.segment_list:
                segment._calculate_mean_dw_mass_per_segment(pipe=pipe1, 
                                        # concentration_drinking_water=concentration_drinking_water_n_min_1,
                                        concentration_drinking_water=pipe1.concentration_drinking_water,
                                        concentration_groundwater=pipe1.concentration_groundwater,)

                sum_mass_segment += segment.mass_chemical_drinkwater

            counter +=1
            criteria = abs(1 - mass_drinkingwater_norm / sum_mass_segment)
            # if abs(1 - mass_drinkingwater_norm / sum_mass_segment) <= tolerance:

            # criteria = abs(1 - concentration_drinking_water_n_min_1 / concentration_drinking_water_n) / relaxation_factor

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

                # two initial guesses to compare the goodness of fit
                if counter == 1:
                    concentration_groundwater_n_plus_1 = pipe1.solubility * 0.9
                if counter == 2:
                    concentration_groundwater_n_plus_1 = pipe1.solubility * 0.5 #['Hexachloorbenzeen', 'indeno[1,2,3-cd]pyrene'] don't work
                if counter >2:
                    if (criteria < criteria_list[counter-1]) or (concentration_groundwater_n_plus_1 < pipe1.concentration_drinking_water):
                        lower_limit = concentration_groundwater_n_min_1
                        concentration_groundwater_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                        proc = 1
                    else:
                        upper_limit = concentration_groundwater_n_min_1
                        concentration_groundwater_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                        proc = 2
                        
                # print(concentration_groundwater_n_min_1, concentration_groundwater_n_plus_1, criteria, proc, lower_limit, upper_limit)
                criteria_old = criteria

            # if counter % 100 ==0 : print(concentration_drinking_water_n_min_1, concentration_drinking_water_n)

new_fails = fails
fails
#%%
# CONCEPT version DW->GW PEAK
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
# database = database.loc[database['log_distribution_coefficient']>=0]
database_chemicals = database['chemical_name_NL']

df = read_csv('did_not_pass_chems_peak_allowable.csv')
failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

# database_chemicals= ['Hexachloorbenzeen']
# chemical_name= 'benzo[b]fluoranthene'

fails = []
over_solubility=[]
for chemical_name in  database_chemicals: #new_fails: #database_chemicals: # new_fails: #failed_chem: # 


    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        # concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    pipe1._fetch_chemical_database(chemical_name=pipe1.chemical_name, 
                                          suppress_print = True, 
                                          language=pipe1.language)

    if pipe1.concentration_drinking_water > pipe1.solubility:
        # raise ValueError(
            over_solubility.append(chemical_name)

            print('Error, the drinking water concentration given or the default drinking water norm is higher than the solubility of the chemical. \
            Input a lower drinking water concentration using .set_conditions()')
    else:
        # calculate initial guess for gw concentration
        sum_KDA_d = 0
        for segment in pipe1.segment_list:
            # calculate the sum of the Kpw * DP * SA *f_stag / d for all pipe segments
            log_Dp_ref = segment._calculate_ref_logD(chemical_group_number=pipe1.chemical_group_number,
                        molecular_weight=pipe1.molecular_weight)
            log_Kpw_ref = segment._calculate_ref_logK(chemical_group_number=pipe1.chemical_group_number,
                        log_octanol_water_partitioning_coefficient=pipe1.log_octanol_water_partitioning_coefficient)
            
            #stagnation factor with reference values for LogDp and LogKpw
            stagnation_factor = 10 ** max((((log_Dp_ref + 12.5) / 2 + 
                                log_Kpw_ref) * 0.73611 + 
                                -1.03574 ), 0)            

            sum_KDA_d_segment = ( 10 ** log_Dp_ref * 10 ** log_Kpw_ref * segment.permeation_surface_area 
                                * stagnation_factor
                                / segment.diffusion_path_length )

            sum_KDA_d += sum_KDA_d_segment

        # initial guess concentration in groundwater
        concentration_groundwater_n_plus_1 = pipe1.concentration_drinking_water * (1 
                                    + pipe1.total_volume * segment.ASSESSMENT_FACTOR_GROUNDWATER 
                                    / pipe1.stagnation_time / sum_KDA_d) 
        
        counter = 0
        max_iterations = 100
        relaxation_factor = 0.5
        tolerance = 0.01
        concentration_drinking_water_n_min_1 = pipe1.concentration_drinking_water
        concentration_drinking_water_n_plus_1 = pipe1.concentration_drinking_water
        lower_limit = pipe1.concentration_drinking_water
        upper_limit = pipe1.solubility
        criteria_list = [0]
        min_criteria = 100
        proc = 0

        while True:  
            concentration_groundwater_n_min_1 = concentration_groundwater_n_plus_1

            pipe1.set_conditions(chemical_name=pipe1.chemical_name,                                    
                concentration_groundwater=concentration_groundwater_n_min_1,
                flow_rate=pipe1.flow_rate,
                concentration_drinking_water=pipe1.concentration_drinking_water, 
                temperature_groundwater=pipe1.temperature_groundwater, 
                stagnation_time = pipe1.stagnation_time,
                suppress_print = True, 
                language = pipe1.language)

            sum_mass_segment = 0

            # mass of chemical in pipe water to meet drinking water norm
            mass_drinkingwater_norm = (pipe1.concentration_drinking_water * pipe1.total_volume)

            for segment in pipe1.segment_list:
                segment._calculate_peak_dw_mass_per_segment(pipe=pipe1, 
                                        # concentration_drinking_water=concentration_drinking_water_n_min_1,
                                        concentration_drinking_water=pipe1.concentration_drinking_water,
                                        concentration_groundwater=pipe1.concentration_groundwater,)

                sum_mass_segment += segment.mass_chemical_drinkwater

            counter +=1
            criteria = abs(1 - mass_drinkingwater_norm / sum_mass_segment)
            # if abs(1 - mass_drinkingwater_norm / sum_mass_segment) <= tolerance:

            # criteria = abs(1 - concentration_drinking_water_n_min_1 / concentration_drinking_water_n) / relaxation_factor

            if criteria <= tolerance:
                print('solution_found!')
                break
            elif counter > max_iterations:
                fails.append(chemical_name)
                print('Max iterations exceeded')
                break
            else:
                # concentration_groundwater_n_plus_1 = (concentration_groundwater 
                #                    * (1 - relaxation_factor + relaxation_factor 
                #                       * (mass_drinkingwater_norm / sum_mass_segment)))
                
                # concentration_groundwater = new_groundwater

                min_criteria = min(min_criteria, criteria)
                criteria_list.append(criteria)

                # two initial guesses to compare the goodness of fit
                if counter == 1:
                    concentration_groundwater_n_plus_1 = pipe1.solubility * 0.99
                if counter == 2:
                    concentration_groundwater_n_plus_1 = pipe1.solubility * 0.1
                if counter >2:
                    if (criteria < criteria_list[counter-1]) or (concentration_groundwater_n_plus_1 < pipe1.concentration_drinking_water):
                        lower_limit = concentration_groundwater_n_min_1
                        concentration_groundwater_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                        proc = 1
                    else:
                        upper_limit = concentration_groundwater_n_min_1
                        concentration_groundwater_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                        proc = 2
                        
                # print(concentration_groundwater_n_min_1, concentration_groundwater_n_plus_1, criteria, proc, lower_limit, upper_limit)
                criteria_old = criteria

            # if counter % 100 ==0 : print(concentration_drinking_water_n_min_1, concentration_drinking_water_n)

len(fails)

new_fails = fails
fails

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
database_chemicals = database['chemical_name_NL']

df = read_csv('did_not_pass_chems_peak_allowable.csv')
failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

for chemical_name in failed_chem: #database_chemicals:
    pipe1.set_conditions(
        chemical_name=chemical_name, 
                        concentration_groundwater =input_gw,
                        temperature_groundwater=12, 
                        flow_rate=0.5)

    pipe1.validate_input_parameters()

    mean_conc=pipe1.calculate_peak_dw_concentration()

    pipe1.set_conditions(chemical_name=chemical_name, 
                        temperature_groundwater=12, 
                        concentration_drinking_water = mean_conc,
                        flow_rate=0.5)

    output_gw = pipe1.calculate_peak_allowable_gw_concentration(relaxation_factor=0.1)

    if abs(1-(input_gw/output_gw)) < 0.01:
        chem_name.append(chemical_name)
        pass_fail.append('passed')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append(output_gw)
    else: 
        chem_name.append(chemical_name)
        pass_fail.append('failed')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append(output_gw)

df_peak = pd.DataFrame(list(zip(chem_name, pass_fail, mean_conc_vals, output_gw_vals)),
               columns =['chemical_name_NL', 'pass_fail', 'mean_conc', 'output_gw'])

database_dict = database.set_index('chemical_name_NL').to_dict()

for key, value in database_dict.items():
    df_peak[key] = df_peak['chemical_name_NL'].map(database_dict[key])

df_peak.to_csv('did_not_pass_chems_peak_allowable_failed.csv')

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
database_chemicals = database['chemical_name_NL']

df = read_csv('did_not_pass_chems_mean_allowable.csv')
failed_chem = list(df['chemical_name_NL'].loc[df.pass_fail == 'failed'])

for chemical_name in database_chemicals:
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

    output_gw = pipe1.calculate_mean_allowable_gw_concentration(relaxation_factor = 0.1)

    if abs(1-(input_gw/output_gw)) < 0.01:
        chem_name.append(chemical_name)
        pass_fail.append('passed')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append(output_gw)
    else: 
        chem_name.append(chemical_name)
        pass_fail.append('failed')
        mean_conc_vals.append(mean_conc)
        output_gw_vals.append(output_gw)

df_mean = pd.DataFrame(list(zip(chem_name, pass_fail, mean_conc_vals, output_gw_vals)),
               columns =['chemical_name_NL', 'pass_fail', 'mean_conc', 'output_gw'])

database_dict = database.set_index('chemical_name_NL').to_dict()

for key, value in database_dict.items():
    df_mean[key] = df_mean['chemical_name_NL'].map(database_dict[key])

df_mean.to_csv('did_not_pass_chems_mean_allowable_failed.csv')

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

