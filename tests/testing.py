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

from pathlib import Path
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory
# This is not working, why??
# try:
#     from project_path import module_path #the dot says look in the current folder, this project_path.py file must be in the folder here
# except ModuleNotFoundError:
#     from project_path import module_path

from ppc.ppc_overview_script import Pipe 

#%%

def test_logKpw_ref():
    '''test against the excel'''

def test_logDp_ref():
    '''test against the excel'''