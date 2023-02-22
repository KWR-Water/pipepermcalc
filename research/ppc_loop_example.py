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

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors

import numpy as np
import pandas as pd
import os
import sys
from pandas import read_csv
from pandas import read_excel
import math
import datetime
from datetime import timedelta
from scipy.optimize import minimize

from project_path import file_path

from pipepermcalc.pipe import * 
from pipepermcalc.segment import * 
# ah_todo need to add an additional import
#%%
# Mean concentration in groundwater for a given drinking water 

seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_flow_rate(flow_rate=0.5)

chemical_name = 'Benzene'
concentration_drinking_water = 0.001 #norm for drinking water
tolerance = 0.01
relaxation_factor = 0.1
max_iterations = 1000
counter = 0

# Calculate the  initial guess for the gw based on the norm/c_dw 

pipe_permeability_dict = pipe1._fetch_chemical_database(
                                chemical_name=chemical_name)

#*** initial guess gw concentration
for segment in pipe1.segment_list:
    log_Dp_ref = segment._calculate_ref_logD(pipe_permeability_dict=pipe_permeability_dict, )
    log_Kpw_ref = segment._calculate_ref_logK(pipe_permeability_dict=pipe_permeability_dict, )

    concentration_groundwater = ((concentration_drinking_water * pipe1.total_volume 
                                  * segment.assessment_factor_groundwater 
                                  * segment.diffusion_path_length) 
                                  / (log_Kpw_ref * log_Dp_ref 
                                     * segment.permeation_surface_area) 
                                     + concentration_drinking_water )

# @martinvdS, which segment to use? or how to sum them?    

while True:
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                temperature_groundwater=12, 
                                concentration_groundwater=concentration_groundwater, 
                                )
    sum_mass_segment = 0

    pipe1.pipe_permeability_dict = pipe1._fetch_chemical_database(chemical_name=pipe1.chemical_name)
    pipe1.pipe_permeability_dict['chemical_name'] = pipe1.chemical_name
    pipe1.pipe_permeability_dict['concentration_groundwater'] = pipe1.concentration_groundwater
    pipe1.pipe_permeability_dict['temperature_groundwater'] = pipe1.temperature_groundwater
    concentration_groundwater = pipe1.pipe_permeability_dict['concentration_groundwater'] 

    for segment in pipe1.segment_list:
        delta_c = concentration_groundwater - concentration_drinking_water
        segment._calculate_pipe_K_D(pipe1.pipe_permeability_dict, 
                                pipe1._groundwater_conditions_set, ) 
        
        segment.mass_chemical_drinkwater = ((segment.permeation_coefficient 
                                             * segment.permeation_surface_area 
                                             * delta_c / segment.diffusion_path_length ) 
                                            / segment.assessment_factor_groundwater)
        

        sum_mass_segment += segment.mass_chemical_drinkwater

    concentration_pipe_drinking_water = (sum_mass_segment / 
                                        pipe1.flow_rate)
    counter +=1

    if abs(1 - concentration_drinking_water / concentration_pipe_drinking_water) <= tolerance:
        break
    elif counter > max_iterations:
        print('Max iterations exceeded')
        break
    else:
        new_groundwater = concentration_groundwater * (1 + ((concentration_drinking_water - concentration_pipe_drinking_water) / concentration_drinking_water * relaxation_factor))
        concentration_groundwater = new_groundwater
        if counter % 100 == 0 : print(new_groundwater)

print(new_groundwater, counter)
#%%
seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

seg2 = Segment(name='seg2',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1 = Pipe(segment_list=[seg1])
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8, 
                                )

pipe1.set_flow_rate(flow_rate=0.5)
pipe1.calculate_mean_dw_concentration()
pipe1.pipe_permeability_dict

#%% COMPLETE
# Mean concentration in drinking water for a given groundwater concentration 

tolerance = 0.1
relaxation_factor = 0.1
max_iterations = 1000

concentration_drinking_water = 0.01 #initial guess
counter = 0

pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8, 
                                )
pipe1.pipe_permeability_dict = pipe1._fetch_chemical_database(chemical_name=pipe1.chemical_name)
pipe1.pipe_permeability_dict['chemical_name'] = pipe1.chemical_name
pipe1.pipe_permeability_dict['concentration_groundwater'] = pipe1.concentration_groundwater
pipe1.pipe_permeability_dict['temperature_groundwater'] = pipe1.temperature_groundwater

concentration_groundwater = pipe1.pipe_permeability_dict['concentration_groundwater'] 

while True:    

    sum_mass_segment = 0

    for segment in pipe1.segment_list:

        # segment function
        delta_c = concentration_groundwater - concentration_drinking_water
        segment._calculate_pipe_K_D(pipe1.pipe_permeability_dict, 
                                pipe1._groundwater_conditions_set, 
                        ) 

        segment.mass_chemical_drinkwater = ((segment.permeation_coefficient * segment.permeation_surface_area * delta_c / segment.diffusion_path_length ) 
                                            / segment.assessment_factor_groundwater)
        # end segment function

        sum_mass_segment += segment.mass_chemical_drinkwater

    concentration_pipe_drinking_water = (sum_mass_segment / 
                                        pipe1.flow_rate) #volume of water consumed in 1 day = flow rate
    counter +=1

    if abs(1 - concentration_drinking_water / concentration_pipe_drinking_water) <= tolerance:
        break
    elif counter > max_iterations:
        print('Max iterations exceeded')
        break
    else:
        concentration_drinking_water = relaxation_factor * concentration_pipe_drinking_water + (1 - relaxation_factor) * concentration_drinking_water

    if counter % 100 ==0 : print(concentration_drinking_water)

print(concentration_pipe_drinking_water, counter)

#%% COMPLETE
# PEAK concentration in drinking water for a given groundwater concentration 
seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1 = Pipe(segment_list=[seg1])
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8, 
                                )

pipe1.set_flow_rate(flow_rate=0.5)

tolerance = 0.1
relaxation_factor = 0.1
max_iterations = 1000

concentration_drinking_water = 0.001 #initial guess
counter = 0

pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8, 
                                )

pipe1.pipe_permeability_dict = pipe1._fetch_chemical_database(chemical_name=pipe1.chemical_name)
pipe1.pipe_permeability_dict['chemical_name'] = pipe1.chemical_name
pipe1.pipe_permeability_dict['concentration_groundwater'] = pipe1.concentration_groundwater
pipe1.pipe_permeability_dict['temperature_groundwater'] = pipe1.temperature_groundwater

concentration_groundwater = pipe1.pipe_permeability_dict['concentration_groundwater'] 
stagnation_time_hours = 8
stagnation_time = stagnation_time_hours / 24

while True:    

    sum_mass_segment = 0

    for segment in pipe1.segment_list:

        # segment function
        delta_c = concentration_groundwater - concentration_drinking_water
        segment._calculate_pipe_K_D(pipe1.pipe_permeability_dict, 
                                pipe1._groundwater_conditions_set, 
                        ) 
        stagnation_factor = segment._calculate_stagnation_factor()

        segment.mass_chemical_drinkwater = ((segment.permeation_coefficient 
                                             * segment.permeation_surface_area 
                                             * delta_c / segment.diffusion_path_length 
                                             * stagnation_time * stagnation_factor) 
                                            / segment.assessment_factor_groundwater)
        # end segment function

        sum_mass_segment += segment.mass_chemical_drinkwater

    concentration_pipe_drinking_water = (sum_mass_segment / 
                                        pipe1.total_volume) 
    counter +=1

    if abs(1 - concentration_drinking_water / concentration_pipe_drinking_water) <= tolerance:
        break
    elif counter > max_iterations:
        print('Max iterations exceeded')
        break
    else:
        concentration_drinking_water = relaxation_factor * concentration_pipe_drinking_water + (1- relaxation_factor) * concentration_drinking_water

    if counter % 100 ==0 : print(concentration_drinking_water)

print(concentration_pipe_drinking_water, counter)

  
#%%
# Peak concentration in groundwater for a given drinking water 

norm_dw = 0.001
tolerance = 0.01
stagnation_time_hours = 8
# How to we prevent that the code goes off in the wrong direction?

initial_guess_gw = 10+

concentration_pipe_drinking_water = 100
new_groundwater = initial_guess_gw
counter = 0

while True:
    guess_gw = new_groundwater
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                temperature_groundwater=12, 
                                concentration_groundwater=guess_gw, #*** initial guess dw concentration
                                )
    sum_mass_segment = 0

    pipe1.pipe_permeability_dict = pipe1._fetch_chemical_database(chemical_name=pipe1.chemical_name)
    pipe1.pipe_permeability_dict['chemical_name'] = pipe1.chemical_name
    pipe1.pipe_permeability_dict['concentration_groundwater'] = pipe1.concentration_groundwater
    pipe1.pipe_permeability_dict['temperature_groundwater'] = pipe1.temperature_groundwater

    for segment in pipe1.segment_list:
        segment._calculate_peak_dw_mass_per_segment(pipe_permeability_dict = pipe1.pipe_permeability_dict,
                                stagnation_time_hours = stagnation_time_hours,  
                                _groundwater_conditions_set = pipe1._groundwater_conditions_set,
                                flow_rate = pipe1.flow_rate                                       
                                )
        sum_mass_segment += segment.mass_chemical_drinkwater
    
    concentration_pipe_drinking_water = (sum_mass_segment / 
                                        pipe1.total_volume)
    counter +=1

    if abs(1 - norm_dw / concentration_pipe_drinking_water) <= tolerance:
        break
    else:
        new_groundwater = initial_guess_gw * (norm_dw / concentration_pipe_drinking_water)

print(new_groundwater, counter)
      








#%%