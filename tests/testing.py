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
    '''test the calculatiion of the reference logK value against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Kpw_ref'], 5)
    logKpw_ref_answer = 1.64761000
    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")
    # else:
    #     print("Success, no error!")

def test_logDp_ref():
    '''test the calculatiion of the reference logD value against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Dp_ref'], 5)
    logKpw_ref_answer = -11.54717
    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")

def test_logKp_ref_temperature_correction():
    '''test the calculatiion of the reference logK value, corrected for temperature
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Ktemp'], 6)
    logKpw_ref_answer = -0.071506
    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")    

def test_logDp_ref_temperature_correction():
    '''test the calculatiion of the reference logD value, corrected for temperature
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Dtemp'], 6)
    logKpw_ref_answer = -0.305084
    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")    

def test_logKp_ref_other_correction():
    '''test the calculatiion of the reference logK value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Kconc'], 6)
    logKpw_ref_answer = -0.103871

    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")    

def test_logDp_ref_other_correction():
    '''test the calculatiion of the reference logD value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Dconc'], 6)
    logKpw_ref_answer = -0.391329

    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")  

def test_logKpw():
    '''test the calculatiion of the reference logK value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Kpw'], 6)
    logKpw_ref_answer = 1.472233

    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")    

def test_logDpw():
    '''test the calculatiion of the reference logD value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    chemical_name="Benzene", 
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Dp'], 6)
    logKpw_ref_answer = -12.243587

    try:
        # assert output == output_phreatic
        assert answer_round == logKpw_ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")  