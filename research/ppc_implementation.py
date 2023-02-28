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
#%%
# # Test the fuzzy-wuzzy matching
# seg1 = Segment(name='seg1',
#                 material='PE40',
#                 length=25,
#                 inner_diameter=0.0196,
#                 thickness=0.0027,
#                 )

# pipe1 = Pipe(segment_list=[seg1])

# database = ["Benzene", "Benzeen", "Methylbenzene", "Toluene", "Ethylbenzene", "Xylene"]

# chemical_name = "Methylbenzeene"

# pipe1._extract_matching_chemical_name(chemical_name=chemical_name, database=database)

#%%

seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1 = Pipe(segment_list=[seg1])
pipe1.set_groundwater_conditions(chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8, 
                                )
pipe1.set_flow_rate(flow_rate=0.5)
# pipe1.calculate_mean_dw_concentration()
pipe1.calculate_peak_dw_concentration()
# pipe1.pipe_permeability_dict

#%%

seg1 = Segment(name='seg1',
            material='PE40',
            length=5,
            inner_diameter=0.0196,
            thickness=0.0027,
            )

seg2 = Segment(name='seg2',
                material='EPDM',
                length=0.06,
                inner_diameter=0.025,
                thickness=0.001,
                diffusion_path_length = 0.06, 
                permeation_direction = 'parallel'
                )

seg3 = Segment(name='seg3',
            material='PE40',
            length=6,
            inner_diameter=0.0196,
            thickness=0.0027,
            )

pipe1 = Pipe(segment_list=[seg1, seg2, seg3])

# pipe1 = Pipe(segment_list=[seg2])
pipe1.set_flow_rate(flow_rate=0.5)

pipe1.set_groundwater_conditions(chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                concentration_groundwater = 1.8)
seg2.__dict__

pipe1.calculate_peak_allowable_gw_concentration(concentration_drinking_water=0.001,
                                chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                tolerance = 0.001
                                )
pipe1.pipe_permeability_dict
#%%
seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1 = Pipe(segment_list=[seg1])
pipe1.segment_list

pipe1.set_groundwater_conditions(chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                concentration_groundwater = 0.112980124482)
pipe1.set_flow_rate(flow_rate=0.5)
pipe1.calculate_peak_dw_concentration()   
pipe1.pipe_permeability_dict


#%%
seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

pipe1 = Pipe(segment_list=[seg1])
pipe1.segment_list

pipe1.set_groundwater_conditions(chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                concentration_groundwater = 0.112980124482)
pipe1.set_flow_rate(flow_rate=0.5)
# pipe1.calculate_peak_dw_concentration()  
# pipe1.pipe_permeability_dict

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
                                chemical_name='Benzeen',)

pipe_permeability_dict['temperature_groundwater'] = 12
pipe_permeability_dict['concentration_groundwater'] = 1.8
pipe_permeability_dict

seg1._calculate_pipe_K_D(pipe_permeability_dict, _groundwater_conditions_set=True)
#%%
pipe1 = Pipe()
pipe1.set_groundwater_conditions(chemical_name="Benzeen", 
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