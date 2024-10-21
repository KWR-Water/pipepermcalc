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
                material='PE40',
                length=25, #meters
                inner_diameter=0.0196, #meters
                wall_thickness=0.0027,) #meters

pipe1 = Pipe(segment_list=[seg1])

pipe1.set_conditions(chemical_name="Benzeen",
                            temperature_groundwater=12, # deg. C
                            # concentration_groundwater=0.1, #g/m3
                            flow_rate=0.5) #m3/day

pipe1.validate_input_parameters()

peak_conc = pipe1.calculate_peak_allowable_gw_concentration()

print("The peak drinking water concentration is:", round(peak_conc,4), "g/m3")

#%%

seg1 = Segment(name='seg1',
            material= 'PE40',
            length=6/1000,
            inner_diameter=23.5/1000,
            wall_thickness=1/1000,
            permeation_direction='parallel',
            diffusion_path_length= 6/1000,
            )

seg2 = Segment(name='seg2',
            material= 'PE40',
            length=0.025 , #/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

seg3 = Segment(name='seg3',
            material= 'PE40',
            length=100/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

seg4 = Segment(name='seg4',
            material= 'PVC',
            length=6,
            inner_diameter=40/1000,
            wall_thickness=2.7/1000,
            )

pipe2 = Pipe(segment_list=[seg1, seg2, seg3, seg4,])

pipe2.set_conditions(chemical_name="Benzeen", 
                            temperature_groundwater=12, # deg. C
                            concentration_groundwater=0.1, #g/m3
                            flow_rate=0.5) #m3/day

pipe2.validate_input_parameters()

mean_conc = pipe2.calculate_mean_dw_concentration()

print("The mean drinking water concentration is:", round(mean_conc,4), "g/m3")

#%%

#%%

seg1 = Segment(name='seg1',
            material= 'PE40',
            length=6/1000,
            inner_diameter=23.5/1000,
            wall_thickness=1/1000,
            permeation_direction='parallel',
            diffusion_path_length= 6/1000,
            )

seg2 = Segment(name='seg2',
            material= 'PE40',
            length=0.025 , #/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

seg3 = Segment(name='seg3',
            material= 'PE40',
            length=100/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

seg4 = Segment(name='seg4',
            material= 'PVC',
            length=6,
            inner_diameter=40/1000,
            wall_thickness=2.7/1000,
            )

pipe2 = Pipe(segment_list=[seg1, seg2, seg3, seg4])
chemicals = ['benzene','ethylbenzene', 'toluene']

for chemical in chemicals:
    with np.errstate(divide='ignore'):
        pipe2.set_conditions(
            concentration_groundwater=0.1, #g/m3
            chemical_name=chemical, 
            temperature_groundwater=12, 
            flow_rate=0.5, 
            suppress_print=True, 
            suppress_warning = True)

    pipe2.validate_input_parameters()

    mean_conc = pipe2.calculate_mean_dw_concentration()

    print("The mean drinking water concentration for", chemical, "is:", round(mean_conc,8), "g/m3")

#%%
#%%

seg1 = Segment(name='seg',
            material= 'PE40',
            length=0.025 , #/1000,
            inner_diameter=33.3/1000,
            wall_thickness=2.7/1000,
            )

pipe3 = Pipe(segment_list=[seg1])
chemicals = ['benzene','ethylbenzene', 'toluene']

for chemical in chemicals:
    pipe3.set_conditions(
            concentration_groundwater=0.1, #g/m3
            chemical_name=chemical, 
            temperature_groundwater=12, 
            flow_rate=0.5, 
            suppress_print=True, 
            suppress_warning = True)

    pipe3.validate_input_parameters()

    mean_conc = pipe3.calculate_mean_dw_concentration()

    print("The mean drinking water concentration for", chemical, "is:", round(mean_conc,8), "g/m3")
#%%
seg1 = Segment(name='seg1', material='PE40', length=25, inner_diameter=0.0196, wall_thickness=0.0027)

pipe3 = Pipe(segment_list=[seg1])

chemicals = ['benzene','ethylbenzene', 'toluene']

for chemical in chemicals:
                 pipe3.set_conditions(
                                 concentration_groundwater=0.1, #g/m3
                                 chemical_name=chemical,
                                 temperature_groundwater=12,
                                 flow_rate=0.5,
                                 suppress_print=True,
                                 suppress_warning = True)
 

pipe3.validate_input_parameters()

mean_conc = pipe3.calculate_mean_dw_concentration()

print("The mean drinking water concentration for", chemical, "is:", round(mean_conc,8), "g/m3")