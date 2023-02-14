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
                diameter=0.0196,
                thickness=0.0027,
                )

# pipe1.add_segment(name='seg2',
#                 material='PE80',
#                 length=25,
#                 diameter=0.0196,
#                 thickness=0.0027,
#                 )


# pipe1.add_segment(name='seg2',
#                 material='PE80',
#                 length=0.001,
#                 diameter=0.0235,
#                 thickness=0.0010,
#                 diffusion_path_length=0.001
#                 )

pipe1.set_flow_rate(flow_rate=0.5)


# pipe1.calculate_peak_dw_concentration(stagnation_time_hours = 8,) 
pipe_segment = 'seg1'
stagnation_time_hours = 8
flow_rate=0.5

stagnation_time = stagnation_time_hours / 24 # days
segment_volume = pipe1.pipe_dictionary['segments'][pipe_segment]['volume']
segment_length = pipe1.pipe_dictionary['segments'][pipe_segment]['length']

segment_surface_area = pipe1.pipe_dictionary['segments'][pipe_segment]['surface_area']
segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 
concentration_groundwater = pipe1.pipe_permeability_dict['concentration_groundwater']
segment_diffusion_path_length = pipe1.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length']
inner_diameter = pipe1.pipe_dictionary['segments'][pipe_segment]['diameter'] - pipe1.pipe_dictionary['segments'][pipe_segment]['thickness'] 
permeation_coefficient = pipe1.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient']

segment_volume = math.pi * (inner_diameter / 2) ** 2 * segment_length

# From equation 4-10 KWR 2016.056
concentratie_drinkwater = ((permeation_coefficient * 2 * concentration_groundwater * 
                           math.pi * (inner_diameter / 2) * segment_length) / 
                            (segment_diffusion_path_length * flow_rate))

#From equation 4-7 in KWR 2016.056
contact_time = (math.pi * (inner_diameter / 2) ** 2 * segment_length) / flow_rate
mass_drinkwater = (permeation_coefficient * concentration_groundwater * 2 * math.pi * (inner_diameter / 2) * segment_length * contact_time / segment_diffusion_path_length)

concentratie_drinkwater2 = mass_drinkwater / segment_volume

concentratie_drinkwater == concentratie_drinkwater2

#%%
concentratie_drinkwater_stagnation = (permeation_coefficient * 2 * concentration_groundwater * stagnation_time) / (segment_diffusion_path_length * inner_diameter/2)


pipe1.pipe_permeability_dict


# %%
