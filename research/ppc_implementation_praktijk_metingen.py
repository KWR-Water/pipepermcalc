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

seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,)

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_conditions(chemical_name="Benzeen",
                            temperature_groundwater=12,
                            concentration_drinking_water=0.1,
                            flow_rate=0.5)

pipe1.validate_input_parameters()

peak_conc = pipe1.calculate_peak_allowable_gw_concentration()

print("The peak groundwater concentration, not exceeding the norm:", round(peak_conc,4), "g/m3")

# mean_conc = pipe1.calculate_mean_dw_concentration()

# print("The mean concentration is:", round(mean_conc,4), "g/m3")
#%%

seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,)

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_conditions(chemical_name="Benzeen",
                            temperature_groundwater=12,
                            concentration_groundwater=peak_conc,
                            flow_rate=0.5)

pipe1.validate_input_parameters()

mean_conc = pipe1.calculate_mean_dw_concentration()

print("The mean concentration is:", round(mean_conc,4), "g/m3")

peak_conc = pipe1.calculate_peak_dw_concentration()

print("The peak concentration is:", round(peak_conc,2), "g/m3")
# %%
