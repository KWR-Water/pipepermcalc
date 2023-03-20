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

mat = 'SBR'

seg1 = Segment(name='seg1',
            material= mat,
            length=7.5/1000,
            inner_diameter=30.3/1000,
            wall_thickness=1.5/1000,
            permeation_direction='parallel',
            diffusion_path_length= 7.5/1000,
            )
seg7 = Segment(name='seg7',
            material= mat,
            length=7.5/1000,
            inner_diameter=30.3/1000,
            wall_thickness=1.5/1000,
            permeation_direction='parallel',
            diffusion_path_length= 7.5/1000,
            )

seg2 = Segment(name='seg2',
            material= mat,
            length=1/1000,
            inner_diameter=23.5/1000,
            wall_thickness=10/1000,
            )

seg3 = Segment(name='seg3',
            material= mat,
            length=6/1000,
            inner_diameter=23.5/1000,
            wall_thickness=1/1000,
            permeation_direction='parallel',
            diffusion_path_length= 6/1000,
            )

seg4 = Segment(name='seg4',
            material= 'PE40',
            length=0.025 , #/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

seg5 = Segment(name='seg5',
            material= 'PE40',
            length=100/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

seg6 = Segment(name='seg6',
            material= 'PVC',
            length=6,
            inner_diameter=40/1000,
            wall_thickness=2.7/1000,
            )

# pipe1 = Pipe(segment_list=[ seg4,])
pipe1 = Pipe(segment_list=[seg1, seg7, seg2, seg3, seg4, seg5, seg6])
# 
chemicals = ['benzene','ethylbenzene', 'toluene']
for chemical in chemicals:
    pipe1.set_conditions(
        chemical_name=chemical, #"fluorene", #
        temperature_groundwater=12, 
        flow_rate=0.5, 
        suppress_print=True )

    pipe1.validate_input_parameters()

    peak_conc = pipe1.calculate_peak_allowable_gw_concentration()
    mean_conc = pipe1.calculate_mean_allowable_gw_concentration()

    print(peak_conc, mean_conc)

#%%

# As the dimensions change, the mean allowable concentration increases, however, 
# the increased concentration is not taken into account in the excel sheet - the 
# log Kpw and Log Dp do not change with the increased concentrations. 
# We can test if this is the issue, by copying the excel supplied mean allowable 
# concentration (C_gws) and inputting that as the gw concentration and calculating the 
# logKpw and LogDp values. These do not match with the excel, and get worse with 
# the decrease in the pipe length. This is because the excel sheet still uses a 
# gw concentration of 1.8 (for benzene) for the calculation of the LogK and LogD, 
# while the code iteratively updates these for each new gw concentration. This 
# leads to differences in the final allowable gw concentration calculated 
# -> excel over estimates the mean allowable gw concentration (see next set of calcs)

# Ultimately the problem is that the PE40 sheets are based on the mean allowable 
# concentration for a 25 m pipe, but the connecting pipes are much, much smaller (1000x).

# The tertiary excel sheet *does* match the PE40 sheet for the mean concentrations, 
# however it does not match for the peak concentrations. The tertiary sheet changes 
# the peak concentration for the changes in the dimensions, while the PE40 sheet 
# does not (which is correct, the dimensions should not change the peak concentration).

# Also the tertiary sheets did not have the assessment factor in them. 
# Also found an error in the peak formula, also a typo elsewhere. fixed now.

Ls = [25, 2.5, 0.25, 0.025]
# mean allowable gw concentrations from the excel sheet PE40 for the lengths above
C_gws = [1.8, 17.988, 179.867, 1798.661]

for L, C_gw in zip (Ls, C_gws):

    # Calculated mean dw conc given the mean allowable gw conc for each length
    # do we get the same dw concentration and Log K and LogD?
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=L , #/1000,
                inner_diameter=19.6/1000,
                wall_thickness=2.7/1000,
                )
    pipe1 = Pipe(segment_list=[ seg1,])
    # C_gw = 1.800

    pipe1.set_conditions(
        concentration_groundwater = C_gw,
        chemical_name='benzeen',
        temperature_groundwater=12, 
        flow_rate=0.5, 
        suppress_print=True )

    pipe1.validate_input_parameters()
    pipe1.calculate_mean_dw_concentration()


    # Same calculation, but fix the logK and logD values to the excel sheet values
    # outcome is that we get the same dw conc as in excel
    concentration_drinking_water =  0.001
    test_Kpw = 1.472233
    test_Dp = -12.243587
    delta_c = C_gw - concentration_drinking_water

    mass_chemical_drinkwater = (((10 ** test_Dp * 10 ** test_Kpw)
                                        * seg1.permeation_surface_area 
                                        * delta_c / seg1.diffusion_path_length ) 
                                        / pipe1.ASSESSMENT_FACTOR_GROUNDWATER
                                        * 24 * 60 * 60)

    mean_dw_conc = mass_chemical_drinkwater / pipe1.flow_rate

    print(mean_dw_conc, 
          pipe1.calculate_mean_dw_concentration(), 
          abs(1-mean_dw_conc/pipe1.calculate_mean_dw_concentration()), 
          test_Dp, 
          seg1.log_Dp)

#%%
# If we calculate the mean allowable gw conc for the different lengths and 
# compare to the excel, we see the excel values are larger than the code values 
# for the very small lengths, therefore the excel over estimates the mean allowable gw concentrations


allowable_gw = []
for L in Ls:

    # Calculated mean allowable gw conc for the different lengths
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=L , #/1000,
                inner_diameter=19.6/1000,
                wall_thickness=2.7/1000,
                )
    pipe1 = Pipe(segment_list=[ seg1,])

    pipe1.set_conditions(
        chemical_name='benzeen',
        temperature_groundwater=12, 
        flow_rate=0.5, 
        suppress_print=True )

    pipe1.validate_input_parameters()
    mean_gw_conc = pipe1.calculate_mean_allowable_gw_concentration()

    allowable_gw.append(mean_gw_conc)
res = [i / j for i, j in zip(allowable_gw, C_gws)]

allowable_gw, C_gws, res


#%%
#MEAN
seg4 = Segment(name='seg4',
            material= 'PE40',
            length=2500/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )
pipe1 = Pipe(segment_list=[ seg4,])

# pipe1 = Pipe(segment_list=[seg1, seg7, seg2, seg3, seg4, seg5, seg6])
pipe1.set_conditions(
    concentration_groundwater=0.161,
    chemical_name='ethylbenzeen',
    temperature_groundwater=12, 
    flow_rate=0.5, 
    suppress_print=True )
pipe1.validate_input_parameters()
pipe1.calculate_peak_dw_concentration()

#%%
mean_mass_drinkingwater = (pipe1.concentration_drinking_water * pipe1.flow_rate)

sum_mass_segment = 0
mass_perc_segment = []
for segment in pipe1.segment_list:
    segment.mass_perc = segment.mass_chemical_drinkwater / mean_mass_drinkingwater
    mass_perc_segment.append(segment.mass_perc*100 )
mass_perc_segment
#%% #PEAK

pipe1 = Pipe(segment_list=[seg1, seg7, seg2, seg3, seg4, seg5, seg6])
pipe1.set_conditions(
    chemical_name='benzeen',
    temperature_groundwater=12, 
    flow_rate=0.5, 
    suppress_print=True )
pipe1.validate_input_parameters()
pipe1.calculate_peak_allowable_gw_concentration()

peak_mass_drinkingwater = (pipe1.concentration_drinking_water * pipe1.total_volume)

sum_mass_segment = 0
mass_perc_segment = []
for segment in pipe1.segment_list:
    segment.mass_perc = segment.mass_chemical_drinkwater / peak_mass_drinkingwater
    mass_perc_segment.append(segment.mass_perc*100 )
mass_perc_segment


#%%

seg2 = Segment(name='seg2',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027)

seg3 = Segment(name='seg3',
                material='PE40',
                length=50,
                inner_diameter=0.0196,
                wall_thickness=0.0027)

pipe2 = Pipe(segment_list=[seg1, seg2])
pipe3 = Pipe(segment_list=[seg3])

pipes = [pipe1, pipe2, pipe3]

for pipe_ in pipes:

    pipe_.set_conditions(concentration_groundwater=1.8,
                        chemical_name="Benzeen", 
                        temperature_groundwater=12,
                        flow_rate=0.5, 
                        suppress_print = True)

    pipe_.validate_input_parameters()

    peak_conc = pipe_.calculate_peak_dw_concentration()
    mean_conc = pipe_.calculate_mean_dw_concentration()

    print(peak_conc, mean_conc)

# Shows that a single segment of PE40 has half as high mean dw concentration as a 
# pipe section which contains 2 segment of the same length. 

# Two segments of the same length gives same mean conc as a single segment which is twice as long. 

# In each case the peak concentration is the same. 

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

        if abs(1-(input_gw/output_gw)) < 0.02:
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

