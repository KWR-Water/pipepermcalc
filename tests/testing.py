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
from project_path import file_path

from pipepermcalc.pipe import * 

#%%

def test_logKpw_ref():
    '''test the calculatiion of the reference logK value against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Kpw_ref'], 5)
    ref_answer = 1.64761000
    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")

def test_logDp_ref():
    '''test the calculatiion of the reference logD value against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Dp_ref'], 5)
    ref_answer = -11.54717
    try:
        
        assert answer_round == ref_answer

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
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Ktemp'], 6)
    ref_answer = -0.071506
    try:
        
        assert answer_round == ref_answer

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
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Dtemp'], 6)
    ref_answer = -0.305084
    try:
        
        assert answer_round == ref_answer

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
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Kconc'], 6)
    ref_answer = -0.103871

    try:
        
        assert answer_round == ref_answer

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
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['f_Dconc'], 6)
    ref_answer = -0.391329

    try:
        
        assert answer_round == ref_answer

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
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Kpw'], 6)
    ref_answer = 1.472233

    try:
        
        assert answer_round == ref_answer

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
                    pipe_material= "PE40",)
    
    answer_round = round(pipe1.pipe_permeability_dict['log_Dp'], 6)
    ref_answer = -12.243587

    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")  

def test_stagnation_factor():
    '''test the calculatiion of the stagnation factor'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    pipe1.set_flow_rate(flow_rate=0.5)
    
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_max_dw_concentration(stagnation_time_hours = 8,pipe_segment = 'seg1')

    answer_round = round(pipe1.pipe_permeability_dict['stagnation_factor'], 6)
    ref_answer = 1.387905

    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.") 

def test_peak_without_stagnation():
    '''test the calculatiion of the peak without stagnation'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    pipe1.set_flow_rate(flow_rate=0.5)
    
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_max_dw_concentration(stagnation_time_hours = 8, 
                                    pipe_segment='seg1',)

    answer_round = round(pipe1.pipe_permeability_dict['concentration_peak_without_stagnation'], 6)
    ref_answer = 0.081403

    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")

def test_peak_with_stagnation():
    '''test the calculatiion of the peak without stagnation'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    pipe1.set_flow_rate(flow_rate=0.5)

    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_max_dw_concentration(stagnation_time_hours = 8, 
                                    pipe_segment='seg1',)

    answer_round = round(pipe1.pipe_permeability_dict['concentration_peak_after_stagnation'], 6)
    ref_answer = 0.112980


    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")

def test_peak_soil_concentration():
    '''test the calculatiion of the peak soil concentration'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    pipe1.set_flow_rate(flow_rate=0.5)

    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_max_dw_concentration(stagnation_time_hours = 8, 
                                    pipe_segment='seg1',)

    answer_round = round(pipe1.pipe_permeability_dict['concentration_peak_soil'], 6)
    ref_answer =   0.171964 

    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")

def test_mean_soil_concentration():
    '''test the calculatiion of the mean soil concentration'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40",)
    
    pipe1.set_flow_rate(flow_rate=0.5)

    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_mean_dw_concentration(
                                    pipe_segment='seg1', 
                                    )
    answer_round = round(pipe1.pipe_permeability_dict['concentration_mean_soil'], 5)
    ref_answer =   2.73921

    try:
        
        assert answer_round == ref_answer

    except AssertionError:
        print("Assertion Exception Raised.")