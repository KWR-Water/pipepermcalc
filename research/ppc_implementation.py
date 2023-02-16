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

# %%
# PLACEHOLDER until the testing function works, ask @Bram
from tests.testing import *

# test1 = test_logKpw_ref()
# test2 = test_logDp_ref()
# test3 = test_logKp_ref_temperature_correction()
# test4 = test_logDp_ref_temperature_correction()
# test5 = test_logKp_ref_other_correction()
# test6 = test_logDp_ref_other_correction()
# test7 = test_logKpw()
# test8 = test_logDpw()
# test9 = test_stagnation_factor()
# test10 = test_peak_without_stagnation()
# test11 = test_peak_with_stagnation()
# test12 = test_peak_soil_concentration()
# test13 = test_mean_soil_concentration()
# test14 = test_updating_partitioning_coefficient()
# test15 = test_updating_diffusion_coefficient()
# test16 = test_calculate_peak_dw_concentration()
# test17 = test_calculate_mean_dw_concentration()
# test18 = test_mean_gw_concentration()
# test19 = test_segment_surface_area_calculations()
#%%
seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

# seg2 = Segment(name='seg2',
#                 material='PE80',
#                 length=25,
#                 inner_diameter=0.0196,
#                 thickness=0.0027,
#                 )

pipe1 = Pipe(segment_list=[seg1])
pipe1.segment_list

pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
pipe1.set_flow_rate(flow_rate=0.5)
pipe1.calculate_mean_dw_concentration()
pipe1.pipe_permeability_dict['mean_concentration_pipe_drinking_water']

#%%
seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

_partitioning_a_dh = 7.92169801506708 #see table 5-6 in KWR 2016.056
_partitioning_b_dh = -17.1875608983359 #see table 5-6 in KWR 2016.056
_diffusion_a_dh = 61.8565740136974 #see table 5-6 in KWR 2016.056
_diffusion_b_dh = -78.9191401984509 #see table 5-6 in KWR 2016.056
assessment_factor_groundwater = 3 
assessment_factor_soil = 1
partitioning_a_c = 0.103965019849463 #see equation 5-20 in KWR 2016.056
partitioning_Cref_Sw = 1.000 #see section 5.4.7 in KWR 2016.056
diffusion_a_c = 0.784077209735583 #see equation 5-18 in KWR 2016.056
diffusion_Cref_Sw = 0.5 #see section 5.4.6 in KWR 2016.056


pipe_permeability_dict = pipe1._fetch_chemical_database(
                                chemical_name='Benzene',)

pipe_permeability_dict['temperature_groundwater'] = 12
pipe_permeability_dict['concentration_groundwater'] = 1.8
pipe_permeability_dict

seg1._calculate_pipe_K_D(pipe_permeability_dict, _groundwater_conditions_set=True)
#%%
pipe1 = Pipe()
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 0.31)
pipe1.add_segment(name='seg1',
                material='PE80',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

# pipe1.add_segment(name='seg2',
#                 material='PE40',
#                 length=25,
#                 inner_diameter=0.0196,
#                 thickness=0.0027,
#                 )

# pipe1.add_segment(name='seg2',
#                 material='PE80',
#                 length=0.001,
#                 inner_diameter=0.0235,
#                 thickness=0.0010,
#                 diffusion_path_length=0.001, )

pipe1.set_flow_rate(flow_rate=0.5)

pipe1.calculate_peak_dw_concentration()
pipe1.pipe_permeability_dict
    

#%%
pipe1 = Pipe()
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 6.858)
# pipe1.add_segment(name='seg1',
#                 material='PE40',
#                 length=25,
#                 inner_diameter=0.0196,
#                 thickness=0.0027,)

pipe1.add_segment(name='seg2',
                material='PE80',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1.set_flow_rate(flow_rate=0.5)

def calculate_peak_allowable_gw_concentration_per_segment(pipe1,drinking_water_norm, stagnation_time_hours,
                                pipe_segment):
    stagnation_time = stagnation_time_hours / 24 # days
    segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']
    segment_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['permeation_surface_area']

    segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 

    #Risk limit value groundwater
    flux_max_stagnation = ( drinking_water_norm * segment_volume /
                    stagnation_time)
    flux_max_stagnation_per_m2 = flux_max_stagnation / segment_surface_area

    stagnation_factor = pipe1._calculate_stagnation_factor(pipe_segment=pipe_segment)

    concentration_gw_peak_without_stagnation = (flux_max_stagnation_per_m2 * 
                                segment_diffusion_path_length / 
                                pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient'] 
                                * pipe1.assessment_factor_groundwater)


    concentration_gw_peak_after_stagnation = stagnation_factor * concentration_gw_peak_without_stagnation
    return concentration_gw_peak_after_stagnation

def calculate_peak_allowable_gw_concentration(drinking_water_norm, stagnation_time_hours = 8, ):

    for pipe_segment in pipe1.pipe_dictionary['segment_list']:
        concentration_gw_peak_after_stagnation = calculate_peak_allowable_gw_concentration_per_segment(pipe1 = pipe1, 
                                                              pipe_segment=pipe_segment,
                                                                     drinking_water_norm=drinking_water_norm,
                                stagnation_time_hours = stagnation_time_hours, )
        print(concentration_gw_peak_after_stagnation)

drinking_water_norm = 0.001 #g/m3

calculate_peak_allowable_gw_concentration(drinking_water_norm=drinking_water_norm, stagnation_time_hours = 8, )

#%%

def objective_function(x, 
                       stagnation_time_hours = 8,
                        pipe_segment = 'seg1',
                        ):
    drinking_water_concentration = x
    drinking_water_norm = pipe1.pipe_permeability_dict['Drinking_water_norm']
    stagnation_time = stagnation_time_hours / 24 # days
    segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']
    segment_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['permeation_surface_area']

    segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 

    stagnation_factor = pipe1._calculate_stagnation_factor(pipe_segment=pipe_segment)

    return (stagnation_factor * (x  * segment_volume /
                stagnation_time / segment_surface_area * 
                            segment_diffusion_path_length / 
                            pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient'] 
                            * pipe1.assessment_factor_groundwater)) - x



# constraints = LinearConstraint(np.ones(n_buyers), lb=n_shares, ub=n_shares)
bounds = [(0,1000)]
res = minimize(objective_function, x0 = 10, bounds = bounds)  #bounds = constraints=constraints)
# X is the initial value of the inputs to the objective function. 
# The objective function is what we’re trying to minimize
answer = res.x[0]

answer

# %%


#%%
def objective_function(x, y=1):
    # y = 1
    return 3 * x ** 4 - 2 * x + y 

# bounds = [(1,3)]
res = minimize(objective_function, 6) #, bounds = bounds)
# X is the initial value of the inputs to the objective function. 
# The objective function is what we’re trying to minimize
answer = res.x[0]
answer
