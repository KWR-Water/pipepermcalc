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

def raise_exception_two_values(answer, ref_answer, round_values=None):
    ''' Raise exception if two values are not equal.'''
    if round_values is None:
        try:
            assert answer == ref_answer
        except AssertionError:
            print("Assertion Exception Raised.")
    else:
        answer_round = round(answer, round_values)
        ref_answer = round(ref_answer, round_values)
        try:
            assert answer_round == ref_answer
        except AssertionError:
            print("Assertion Exception Raised.")


def test_logKpw_ref():
    '''test the calculatiion of the reference logK value against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['log_Kpw_ref'], 
                               ref_answer = 1.64761000, 
                               round_values=5)

def test_logDp_ref():
    '''test the calculatiion of the reference logD value against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
   
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['log_Dp_ref'], 
                               ref_answer = -11.54717, 
                               round_values=5)

def test_logKp_ref_temperature_correction():
    '''test the calculatiion of the reference logK value, corrected for temperature
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['f_Ktemp'], 
                               ref_answer = -0.071506, 
                               round_values=6)


def test_logDp_ref_temperature_correction():
    '''test the calculatiion of the reference logD value, corrected for temperature
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['f_Dtemp'], 
                               ref_answer = -0.305084,
                               round_values=6)
    

def test_logKp_ref_other_correction():
    '''test the calculatiion of the reference logK value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['f_Kconc'],
                               ref_answer = -0.103871,
                               round_values=6)
    

def test_logDp_ref_other_correction():
    '''test the calculatiion of the reference logD value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['f_Dconc'],
                               ref_answer =  -0.391329, 
                               round_values=6)
    

def test_logKpw():
    '''test the calculatiion of the reference logK value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['log_Kpw'],
                               ref_answer = 1.472233,
                               round_values=6)
    

def test_logDpw():
    '''test the calculatiion of the reference logD value, 
    corrected for ?? @MartinvdS column AC
      against the excel'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['log_Dp'], 
                               ref_answer = -12.243587, 
                               round_values=6)
    

def test_stagnation_factor():
    '''test the calculatiion of the stagnation factor'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
   
    pipe1.set_flow_rate(flow_rate=0.5)
    
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_peak_allowable_gw_concentration(stagnation_time_hours = 8)

    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['stagnation_factor'],
                               ref_answer =  1.387905, 
                               round_values=6)


def test_peak_without_stagnation():
    '''test the calculatiion of the peak without stagnation'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    
    pipe1.set_flow_rate(flow_rate=0.5)
    
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_peak_allowable_gw_concentration(stagnation_time_hours = 8, )

    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['concentration_peak_without_stagnation'], 
                               ref_answer = 0.081403, 
                               round_values=6)

def test_peak_with_stagnation():
    '''test the calculatiion of the peak without stagnation'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    
    pipe1.set_flow_rate(flow_rate=0.5)

    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_peak_allowable_gw_concentration(stagnation_time_hours = 8,)
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['concentration_peak_after_stagnation'], 
                               ref_answer = 0.112980, 
                               round_values=6)


def test_peak_soil_concentration():
    '''test the calculatiion of the peak soil concentration'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    
    pipe1.set_flow_rate(flow_rate=0.5)

    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_peak_allowable_gw_concentration(stagnation_time_hours = 8,)
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['concentration_peak_soil'],
                               ref_answer =  0.171964 , 
                               round_values=6)


def test_mean_soil_concentration():
    '''test the calculatiion of the mean soil concentration'''

    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    
    pipe1.set_flow_rate(flow_rate=0.5)

    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1.calculate_mean_allowable_gw_concentration()
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['concentration_mean_soil'], 
                               ref_answer = 2.73921, 
                               round_values=5)


def test_updating_partitioning_coefficient():
    ''' Test the update function for the partitioning coefficient '''
    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1._update_partitioning_coefficient(new_log_Kpw= 0.9116730996845103, 
                                        segment_name='seg1')
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['log_Kpw'],
                               ref_answer = 0.911673, 
                               round_values=6)

    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['permeation_coefficient'], 
                               ref_answer = 4.023463562623052e-07, 
                               round_values=None)


def test_updating_diffusion_coefficient():
    ''' Test the update function for the diffusion coefficient '''
    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )

    pipe1._update_diffusion_coefficient(new_log_Dp= -12.743586769549616, 
                                        segment_name='seg1')
    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['log_Dp'], 
                               ref_answer = -12.743586769549616, 
                               round_values=None)

    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['segments']['seg1']['permeation_coefficient'], 
                               ref_answer = 4.6255147758415636e-07, 
                               round_values=None)


def test_calculate_peak_dw_concentration():
    ''' Test the calculation for the peak concentration in drinking water given 
    a groundwater concentration '''
    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    pipe1.set_flow_rate(flow_rate=0.5)
    pipe1.calculate_peak_dw_concentration()

    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['peak_concentration_pipe_drinking_water'], 
                               ref_answer = 0.001, 
                               round_values=5)

def test_calculate_mean_dw_concentration():
    ''' Test the calculation for the mean concentration in drinking water given 
    a groundwater concentration '''
    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8)
    pipe1.add_segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
    pipe1.set_flow_rate(flow_rate=0.5)
    pipe1.calculate_mean_dw_concentration()

    raise_exception_two_values(answer=pipe1.pipe_permeability_dict['mean_concentration_pipe_drinking_water'], 
                               ref_answer = 0.001, 
                               round_values=5)


