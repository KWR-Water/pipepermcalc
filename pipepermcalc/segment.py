#%% ----------------------------------------------------------------------------
# A. Hockin, January 2023
# KWR 403230-003
# Pipe permeation calculator
# With Martin vd Schans, Bram Hillebrand, Lennart Brokx
#
# ------------------------------------------------------------------------------
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

class segment(name=None,
            material=None,
            length=None,
            inner_diameter=None,
            thickness=None,
            permeation_direction='perpendicular',
            diffusion_path_length=None, ):
    ''' 
    Segment object class to make segments of the pipe 

    Attributes
    -------
    '''

    count = 0 # count of pipe segments

    def __init__(self, 
                ):
        '''
        
        '''  

    def add_segment():
        ''''''       
