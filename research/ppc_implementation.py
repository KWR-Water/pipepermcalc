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
#                 material='PE80',
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

pipe1.calculate_peak_dw_concentration(stagnation_time_hours = 8,)
concentration_drinkwater = pipe1.pipe_permeability_dict['peak_concentration_pipe_drinking_water']
#%%
pipe1.calculate_mean_dw_concentration()
concentration_drinkwater = pipe1.pipe_permeability_dict['mean_concentration_pipe_drinking_water']

round(concentration_drinkwater, 6)
#%% 
pipe_segment = 'seg1'
stagnation_time_hours = 8
flow_rate=0.5


stagnation_time = stagnation_time_hours / 24 # days
segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']

segment_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['inner_surface_area']
segment_length = pipe1.pipe_dictionary['segments'][pipe_segment]['length']
segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 
concentration_groundwater = pipe1.pipe_permeability_dict['concentration_groundwater']
segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length']
inner_diameter = pipe1.pipe_dictionary['segments'][pipe_segment]['inner_diameter'] 
permeation_coefficient = pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient']
flow_rate = pipe1.flow_rate

# @MartinvdS how are the assessment factor incorporated?
# these are different formulas than used in the excel sheet...?

# From equation 4-10 KWR 2016.056
concentration_drinkwater = ((permeation_coefficient * 2 * concentration_groundwater * stagnation_time) / 
                    (segment_diffusion_path_length * (inner_diameter / 2)))

# Add the assessment factor and stagnation factor here @ah_todo
stagnation_factor = pipe1._calculate_stagnation_factor(pipe_segment=pipe_segment)
assessment_factor = pipe1.assessment_factor_groundwater

concentration_drinkwater = concentration_drinkwater / stagnation_factor / assessment_factor

#From equation 4-7 in KWR 2016.056
contact_time = stagnation_time
mass_drinkwater = concentration_drinkwater * segment_volume 

# mass_drinkwater = mass_drinkwater / stagnation_factor / assessment_factor

pipe1.pipe_permeability_dict['segments'][pipe_segment]['concentration_drinkwater'] = concentration_drinkwater
pipe1.pipe_permeability_dict['segments'][pipe_segment]['contact_time'] = contact_time
pipe1.pipe_permeability_dict['segments'][pipe_segment]['mass_drinkwater'] = mass_drinkwater
pipe1.pipe_permeability_dict['segments'][pipe_segment]['stagnation_factor'] = stagnation_factor

round(concentration_drinkwater, 5)
#%%

# C_dw = ( permeation_coefficient * concentration_groundwater * 2  * math.pi * (inner_diameter / 2) * segment_length ) / (segment_diffusion_path_length * flow_rate)

C_dw = ( permeation_coefficient * concentration_groundwater * math.pi * inner_diameter * segment_length) / (segment_diffusion_path_length * flow_rate * assessment_factor)

# C_dw = 0.001

# C_gw = (segment_diffusion_path_length * stagnation_factor * assessment_factor) / (concentration_groundwater * permeation_coefficient * 2 * stagnation_time) 


flux_j = (C_dw * flow_rate)
flux_J = C_dw * flow_rate / segment_surface_area

# C_gw = flux_J * segment_diffusion_path_length * assessment_factor / permeation_coefficient + C_dw
# C_gw = C_dw * flow_rate / segment_surface_area * segment_diffusion_path_length * assessment_factor / permeation_coefficient + C_dw

C_gw = 1.8

# let Cdw be szero
C_dw = ( permeation_coefficient  * C_gw * segment_surface_area) / (flow_rate * segment_diffusion_path_length * assessment_factor) 

flux_J, flux_j, C_dw
# C_dw

# %%
concentration_drinkwater = ((permeation_coefficient *  concentration_groundwater * 
                    segment_surface_area) / 
                    (segment_diffusion_path_length * flow_rate * pipe1.assessment_factor_groundwater ))


#From equation 4-7 in KWR 2016.056
contact_time = (math.pi * (inner_diameter / 2) ** 2 * segment_length) / flow_rate
mass_drinkwater = (permeation_coefficient * concentration_groundwater * 2 * math.pi * (inner_diameter / 2) * segment_length * contact_time / segment_diffusion_path_length)

mass_drinkwater / segment_volume / assessment_factor