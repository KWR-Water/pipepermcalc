#%% ----------------------------------------------------------------------------
# A. Hockin, January 2023
# KWR 403230-003
# Pipe permeation calculator
# With Martin vd Schans, Bram
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

from pipepermcalc.ppc_overview_script import * 

pipe1 = Pipe()
pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
pipe1.calculate_pipe_K_D(
                   chemical_name="Benzene", 
                   pipe_material= "PE40",
                    )
pipe1.pipe_permeability_dict


# pipe1 = Pipe()
# pipe1.add_segment(name='seg1',
#                 material='PE40',
#                 length=25,
#                 diameter=0.0196,
#                 thickness=0.0027,
#                 )
# dict1 = pipe1.pipe_dictionary

# pipe1.add_segment(name='seg2',
#                 material='PE80',
#                 length=25,
#                 diameter=1,
#                 thickness=0.03,
#                 )
# # dict2 = pipe1.pipe_dictionary

# pipe1.pipe_dictionary
#%%

