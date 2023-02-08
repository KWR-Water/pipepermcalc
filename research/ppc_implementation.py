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
pipe1.add_segment(name='seg2',
                material='PE80',
                length=2,
                diameter=0.0196,
                thickness=0.0027,
                )

pipe1.set_flow_rate(flow_rate=0.5)

pipe1.calculate_pipe_K_D(
                   pipe_material= "PE40",
                    )
# pipe1.pipe_dict['seg1']['K'] = 21 # ah_todo add this as example to read the docs in advanced user tutorial
# print(pipe1)
pipe1.pipe_permeability_dict

#%%


# pipe1.pipe_dictionary

pipe1.calculate_max_dw_concentration(stagnation_time_hours = 8, 
                                    pipe_segment='seg1', 
                                    )
pipe1.calculate_mean_dw_concentration(
                                    pipe_segment='seg1', 
                                    )

pipe1.pipe_permeability_dict

#%%
pipe1.calculate_pipe_K_D_per_segment()
pipe1.pipe_permeability_dict

# %%
