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
# ah_todo need to add an additional import

# %%
# PLACEHOLDER until the testing function works, ask @Bram
from tests.testing import *

test1 = test_logKpw_ref()
test2 = test_logDp_ref()
test3 = test_logKp_ref_temperature_correction()
test4 = test_logDp_ref_temperature_correction()
test5 = test_logKp_ref_other_correction()
test6 = test_logDp_ref_other_correction()
test7 = test_logKpw()
test8 = test_logDpw()
test9 = test_stagnation_factor()
test10 = test_peak_without_stagnation()
test11 = test_peak_with_stagnation()
test12 = test_peak_soil_concentration()
test13 = test_mean_soil_concentration()
test14 = test_updating_partitioning_coefficient()
test15 = test_updating_diffusion_coefficient()
test16 = test_calculate_peak_dw_concentration()
test17 = test_calculate_mean_dw_concentration()
test18 = test_mean_gw_concentration()
test19 = test_segment_surface_area_calculations()
#%%
pipe1 = Pipe()
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
pipe1.add_segment(name='seg1',
                material='PE40',
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
#                 diffusion_path_length=0.001
#                 )
pipe1.set_flow_rate(flow_rate=0.5)

pipe1.calculate_mean_dw_concentration()
pipe1.pipe_permeability_dict

#%%

#%%
pipe1 = Pipe()
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
pipe1.add_segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1.set_flow_rate(flow_rate=0.5)

# pipe1.calculate_mean_dw_concentration()
pipe1.pipe_permeability_dict
pipe_segment = 'seg1'

segment_length = pipe1.pipe_dictionary['segments'][pipe_segment]['length']
# segment_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['inner_surface_area']
segment_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['permeation_surface_area']
segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']

segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 
concentration_groundwater = pipe1.pipe_permeability_dict['concentration_groundwater']
segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length']
segment_inner_diameter = pipe1.pipe_dictionary['segments'][pipe_segment]['inner_diameter'] 
permeation_coefficient = pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient']
flow_rate = pipe1.flow_rate

#Excel
concentration_drinkwater_excel = (concentration_groundwater * (segment_surface_area * permeation_coefficient) / 
                                    ((flow_rate * segment_diffusion_path_length * pipe1.assessment_factor_groundwater) + segment_surface_area * permeation_coefficient ))
mass_drinkwater_excel = concentration_drinkwater_excel * segment_volume
mass_flux_excel = (concentration_drinkwater_excel * flow_rate) / segment_surface_area
mass_drinkwater_excel

#Report
concentration_drinkwater = 0.001
contact_time = (math.pi * (segment_inner_diameter / 2) ** 2 * segment_length) / flow_rate
contact_time2 = (segment_volume) / flow_rate
mass_flux = permeation_coefficient * (concentration_groundwater - concentration_drinkwater) / (segment_diffusion_path_length * pipe1.assessment_factor_groundwater)

mass_drinkwater_report = mass_flux * segment_surface_area * contact_time2
mass_drinkwater_excel/mass_drinkwater_report

# #calculation rearranged
mass_drinkwater_calc = (concentration_groundwater * segment_volume * segment_surface_area * permeation_coefficient) / (segment_diffusion_path_length * pipe1.assessment_factor_groundwater * flow_rate + permeation_coefficient * segment_surface_area)
mass_drinkwater_calc
concentration_drinkwater_calc = mass_drinkwater_calc / segment_volume
concentration_drinkwater_calc
# mass_flux, mass_flux_excel #Mass fluxes ok.
#%%

def objective_function(x, 
                       stagnation_time_hours = 8,
                        pipe_segment = 'seg1',
                        ):
    drinking_water_concentration = x
    stagnation_time = stagnation_time_hours / 24 # days
    segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']
    segment_inner_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['inner_surface_area']
    segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 

    flux_max_stagnation = ( drinking_water_concentration * segment_volume /
                stagnation_time)
    flux_max_stagnation_per_m2 = flux_max_stagnation / segment_inner_surface_area

    stagnation_factor = pipe1._calculate_stagnation_factor(pipe_segment=pipe_segment)

    concentration_gw_peak_without_stagnation = (flux_max_stagnation_per_m2 * 
                                segment_diffusion_path_length / 
                                pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient'] 
                                * pipe1.assessment_factor_groundwater)


    concentration_gw_peak_after_stagnation = stagnation_factor * concentration_gw_peak_without_stagnation    
    return concentration_gw_peak_after_stagnation

# constraints = LinearConstraint(np.ones(n_buyers), lb=n_shares, ub=n_shares)
res = minimize(objective_function, x0 = 6)  #bounds = constraints=constraints)
# X is the initial value of the inputs to the objective function. 
# The objective function is what we’re trying to minimize
answer = res.x[0]
answer

#%%
stagnation_time_hours = 8
pipe_segment = 'seg1'
drinking_water_norm = pipe1.pipe_permeability_dict['Drinking_water_norm']

drinking_water_concentration = drinking_water_norm / 1000
stagnation_time = stagnation_time_hours / 24 # days
segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']
segment_inner_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['inner_surface_area']
segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 

#Risk limit value groundwater
flux_max_stagnation = ( drinking_water_concentration * segment_volume /
                stagnation_time)
flux_max_stagnation_per_m2 = flux_max_stagnation / segment_inner_surface_area

stagnation_factor = pipe1._calculate_stagnation_factor(pipe_segment=pipe_segment)

concentration_gw_peak_without_stagnation = (flux_max_stagnation_per_m2 * 
                            segment_diffusion_path_length / 
                            pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient'] 
                            * pipe1.assessment_factor_groundwater)


concentration_gw_peak_after_stagnation = stagnation_factor * concentration_gw_peak_without_stagnation

#Risk limit value soil, first check if a distribution coefficient is known
if math.isnan(pipe1.pipe_permeability_dict['log_distribution_coefficient']):
    concentration_peak_soil = 'log_distribution_coefficient (Kd) unknown'
else:
    concentration_peak_soil = (10 ** pipe1.pipe_permeability_dict['log_distribution_coefficient'] * 
                                concentration_gw_peak_after_stagnation * 
                                pipe1.assessment_factor_soil / pipe1.assessment_factor_groundwater)

concentration_gw_peak_after_stagnation
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['stagnation_time_hours'] = stagnation_time_hours
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['flux_max_stagnation'] = flux_max_stagnation
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['flux_max_stagnation_per_m2'] = flux_max_stagnation_per_m2
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['stagnation_factor'] = stagnation_factor
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['concentration_gw_peak_without_stagnation'] = concentration_gw_peak_without_stagnation
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['concentration_gw_peak_after_stagnation'] = concentration_gw_peak_after_stagnation
# pipe1.pipe_permeability_dict['segments'][pipe_segment]['concentration_peak_soil'] = concentration_peak_soil

#%%


#%%
def objective_function(x):

    return 3 * x ** 4 - 2 * x + 1

res = minimize(objective_function, 6)
# X is the initial value of the inputs to the objective function. 
# The objective function is what we’re trying to minimize
answer = res.x[0]
answer
