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

from pipepermcalc.pipe import * 
from pipepermcalc.segment import * 

#%%

def raise_exception_two_values(answer, ref_answer, round_values=None):
    ''' Raise exception if two values are not equal.'''
    if round_values is None:
        assert answer == ref_answer
    else:
        answer_round = round(answer, round_values)
        ref_answer = round(ref_answer, round_values)
        assert answer_round == ref_answer


def test_logKpw_ref():
    '''test the calculatiion of the reference logK value against the excel'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482, 
                                    flow_rate=0.5)
    raise_exception_two_values(answer=seg1.log_Kpw_ref, 
                               ref_answer = 1.64761000, 
                               round_values=5)

def test_logDp_ref():
    '''test the calculatiion of the reference logD value against the excel'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482, 
                                    flow_rate=0.5)
  
    raise_exception_two_values(answer=seg1.log_Dp_ref, 
                               ref_answer = -11.54717, 
                               round_values=5)

def test_logKp_ref_temperature_correction():
    '''test the calculatiion of the reference logK value, corrected for temperature
      against the excel'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482, 
                                    flow_rate=0.5)
    raise_exception_two_values(answer=seg1.f_Ktemp, 
                               ref_answer = -0.071506, 
                               round_values=6)


def test_logDp_ref_temperature_correction():
    '''test the calculatiion of the reference logD value, corrected for temperature
      against the excel'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482, 
                                    flow_rate=0.5)
    raise_exception_two_values(answer=seg1.f_Dtemp, 
                               ref_answer = -0.305084,
                               round_values=6)
    

def test_logKp_ref_concentration_correction():
    '''test the calculation of the reference logK value, 
    corrected for concentration against the excel'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    
    raise_exception_two_values(answer=seg1.f_Kconc,
                               ref_answer = -0.103871,
                               round_values=6)
    

def test_logDp_ref_concentration_correction():
    '''test the calculatiion of the reference logD value, 
    corrected for concentration '''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater =1.8,
                                    flow_rate=0.5)
    raise_exception_two_values(answer=seg1.f_Dconc,
                               ref_answer =  -0.391329, 
                               round_values=6)
    

def test_logKpw():
    '''test the calculatiion of the logKpw'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    raise_exception_two_values(answer=seg1.log_Kpw,
                               ref_answer = 1.472233,
                               round_values=6)
    

def test_logDpw():
    '''test the calculatiion of the logDw'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    raise_exception_two_values(answer=seg1.log_Dp, 
                               ref_answer = -12.243587, 
                               round_values=6)
    

def test_stagnation_factor():
    '''test the calculatiion of the stagnation factor'''

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    
    pipe1.validate_input_parameters()
    pipe1.calculate_peak_dw_concentration()    

    raise_exception_two_values(answer=seg1.stagnation_factor,
                               ref_answer =  1.387905, 
                               round_values=6)

def test_updating_partitioning_coefficient():
    ''' Test the update function for the partitioning coefficient '''
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    pipe1.validate_input_parameters()
    seg1.log_Kpw= 0.9116730996845103

    pipe1.calculate_peak_dw_concentration()    

    raise_exception_two_values(answer=seg1.log_Kpw,
                               ref_answer = 0.911673, 
                               round_values=6)


def test_updating_diffusion_coefficient():
    ''' Test the update function for the diffusion coefficient '''
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    pipe1.validate_input_parameters()
    seg1.log_Dp= -12.743586769549616
    
    pipe1.calculate_peak_dw_concentration()    
   
    raise_exception_two_values(answer=seg1.log_Dp, 
                               ref_answer = -12.743586769549616, 
                               round_values=None)



def test_calculate_peak_dw_concentration():
    ''' Test the calculation for the peak concentration in drinking water given 
    a groundwater concentration '''
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482,
                                    flow_rate=0.5)
    pipe1.validate_input_parameters()
    
    pipe1.calculate_peak_dw_concentration()    

    raise_exception_two_values(answer=pipe1.concentration_drinking_water, 
                               ref_answer = 0.001, 
                               round_values=3)
    
def test_calculate_peak_dw_mass():
    ''' Test the calculation for the peak concentration in drinking water given 
    a groundwater concentration '''
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 0.112980124482,
                                    flow_rate=0.5)
    pipe1.validate_input_parameters()
    
    pipe1.calculate_peak_dw_concentration()    

    raise_exception_two_values(answer=seg1.mass_chemical_drinkwater, 
                               ref_answer = 7.535e-06, 
                               round_values=9)    

def test_calculate_mean_dw_concentration():
    ''' Test the calculation for the mean concentration in drinking water given 
    a groundwater concentration '''
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    pipe1.validate_input_parameters()
    
    pipe1.calculate_mean_dw_concentration()    

    raise_exception_two_values(answer=pipe1.concentration_drinking_water, 
                               ref_answer = 0.001, 
                               round_values=5)

def test_calculate_mean_dw_mass():
    ''' Test the calculation for the mean concentration in drinking water given 
    a groundwater concentration '''
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    wall_thickness=0.0027,
                    )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list

    pipe1.set_conditions(chemical_name="Benzeen", 
                                    temperature_groundwater=12, 
                                    concentration_groundwater = 1.8,
                                    flow_rate=0.5)
    pipe1.validate_input_parameters()
    
    pipe1.calculate_mean_dw_concentration()    

    raise_exception_two_values(answer=seg1.mass_chemical_drinkwater, 
                               ref_answer = 0.0005000956341327644, 
                               round_values=5)    

def test_segment_surface_area_calculations():
    ''' Test the calculation for the different surface area options '''
    seg1 = Segment(name='seg1',
                material='PE40',
                length=7.5/1000,
                inner_diameter=30.3/1000,
                wall_thickness=1.5/1000,
                permeation_direction='parallel',
                diffusion_path_length=7.5/1000,
                    )

   
    raise_exception_two_values(answer=seg1.permeation_surface_area, 
                               ref_answer = 0.000073159839, 
                               round_values=12)
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=1/1000,
                    inner_diameter=28.5/1000,
                    wall_thickness=10/1000,
                    permeation_direction='perpendicular',
                    diffusion_path_length=10/1000
                    )

    
    raise_exception_two_values(answer=seg1.permeation_surface_area, 
                               ref_answer = 0.000089535391, 
                               round_values=12)

    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=33.3/1000,
                    inner_diameter=25/1000,
                    wall_thickness=2.7/1000,
                    permeation_direction='perpendicular',)
    
    raise_exception_two_values(answer=seg1.permeation_surface_area, 
                               ref_answer = 0.002615375884, 
                               round_values=12)


def test_return_matching_chemical_name():
    seg1 = Segment(name='seg1',
                material='PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.segment_list
    chemical_name = 'Benzeen'
    matching_name = pipe1._extract_matching_chemical_name(chemical_name=chemical_name, 
                                          database=list(pipe1.ppc_database['chemical_name_NL']))
    assert chemical_name == matching_name


def test_calculate_mean_allowable_gw_concentration():
    ''' Test the calculation for the mean concentration allowed in groundwater given 
    a drinking water concentration '''
    
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.set_conditions(
        chemical_name="Benzeen", 
        temperature_groundwater=12, 
        concentration_drinking_water=0.001,
        flow_rate=0.5 )

    pipe1.validate_input_parameters()

    mean_conc = pipe1.calculate_mean_allowable_gw_concentration(tolerance = 0.001)

    raise_exception_two_values(answer=mean_conc, 
                               ref_answer = 1.8011, 
                               round_values=4)

def test_calculate_peak_allowable_gw_concentration():
    ''' Test the calculation for the peak concentration allowed in groundwater given 
    a drinking water concentration '''
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])
    pipe1.set_conditions(
        chemical_name="Benzeen", 
        temperature_groundwater=12, 
        concentration_drinking_water=0.001,
        flow_rate=0.5 )

    pipe1.validate_input_parameters()

    peak_conc = pipe1.calculate_peak_allowable_gw_concentration(tolerance = 0.01)

    raise_exception_two_values(answer=peak_conc, 
                               ref_answer = 0.1124, 
                               round_values=4)

def test_groundwater_to_soil_conversion():
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(
        chemical_name="Benzeen", 
        temperature_groundwater=12, 
        concentration_groundwater=1.8,
        flow_rate=0.5 )

    pipe1.validate_input_parameters()

    pipe1.calculate_mean_dw_concentration()

    raise_exception_two_values(answer=pipe1.concentration_soil, 
                               ref_answer =2.74, 
                               round_values=3)
    # Check that pipe summation is going correctly
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(
        chemical_name="ethylbenzene", 
        temperature_groundwater=12, 
        concentration_groundwater=0.277,
        flow_rate=0.5 )

    pipe1.validate_input_parameters()

    pipe1.calculate_mean_dw_concentration()

    raise_exception_two_values(answer=pipe1.concentration_soil, 
                               ref_answer =1.842, 
                               round_values=3)

def test_soil_to_groundwater_conversion():
    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])

    pipe1.set_conditions(
        chemical_name="Benzeen",
        temperature_groundwater=12, 
        concentration_soil=2.7527429729399238,
        flow_rate=0.5 )

    pipe1.validate_input_parameters()

    pipe1.calculate_mean_dw_concentration()

    raise_exception_two_values(answer=pipe1.concentration_groundwater, 
                               ref_answer = 1.8085521338873323, 
                               round_values=5)

def test_GW_to_DW_to_GW_peak():
    '''Tests if the calculation from GW to DW gives the same result as DW to GW 
    for the peak concentrations'''

    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])
    input_gw = 1

    database = pipe1.view_database_chemical_names( language='NL')
    database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])
    database = database.loc[database['log_distribution_coefficient']>=0]
    database = database.loc[database['Drinking_water_norm'] < database['solubility'] ]
    database_chemicals = database['chemical_name_NL']
    solubilities = database['solubility']

    failed = []

    for chemical_name, solubiliy in zip(database_chemicals, solubilities):
        if input_gw > solubiliy:
            input_gw = 0.01 * solubiliy

        pipe1.set_conditions(
            chemical_name=chemical_name, 
                            concentration_groundwater =input_gw,
                            temperature_groundwater=12, 
                            flow_rate=0.5)

        pipe1.validate_input_parameters()

        peak_conc=pipe1.calculate_peak_dw_concentration()

        pipe1.set_conditions(chemical_name=chemical_name, 
                            temperature_groundwater=12, 
                            concentration_drinking_water = peak_conc,
                            flow_rate=0.5)

        output_gw = pipe1.calculate_peak_allowable_gw_concentration()

        if abs(1-(input_gw/output_gw)) < 0.02:
            pass
        else: 
            failed.append(chemical_name)
    
    assert len(failed) == 0


def test_GW_to_DW_to_GW_mean():
    '''Tests if the calculation from GW to DW gives the same result as DW to GW 
    for the mean concentrations'''

    seg1 = Segment(name='seg1',
                material= 'PE40',
                length=25,
                inner_diameter=0.0196,
                wall_thickness=0.0027,
                )

    pipe1 = Pipe(segment_list=[seg1])
    input_gw = 1

    database = pipe1.view_database_chemical_names( language='NL')
    database = pipe1.ppc_database.dropna(subset=['molecular_weight', 'solubility', 'Drinking_water_norm'])
    database = database.loc[database['log_distribution_coefficient']>=0]
    database = database.loc[database['Drinking_water_norm'] < database['solubility'] ]
    database_chemicals = database['chemical_name_NL']
    solubilities = database['solubility']

    failed = []

    for chemical_name, solubiliy in zip(database_chemicals, solubilities):
        if input_gw > solubiliy:
            input_gw = 0.01 * solubiliy

        pipe1.set_conditions(
            chemical_name=chemical_name, 
                            concentration_groundwater =input_gw,
                            temperature_groundwater=12, 
                            flow_rate=0.5)

        pipe1.validate_input_parameters()

        mean_conc=pipe1.calculate_mean_dw_concentration()

        pipe1.set_conditions(chemical_name=chemical_name, 
                            temperature_groundwater=12, 
                            concentration_drinking_water = mean_conc,
                            flow_rate=0.5)

        output_gw = pipe1.calculate_mean_allowable_gw_concentration()

        if abs(1-(input_gw/output_gw)) < 0.02:
            pass
        else: 
            failed.append(chemical_name)
    
    assert len(failed) == 0

#%%
# Pipe functions
#  _validate_object
# validate_input_parameters
#? _fuzzy_min_score
#* _extract_matching_chemical_name
#? set_conditions
#? _fetch_chemical_database
#* calculate_mean_dw_concentration
#* calculate_peak_dw_concentration
#* calculate_mean_allowable_gw_concentration
#* calculate_peak_allowable_gw_concentration

#Segment functions
#* _correct_for_temperature
#* _concentration_correction
# na _correct_for_age
#* _calculate_ref_logK
#* _calculate_ref_logD
#* _calculate_logK
#* _calculate_logD
#* _calculate_pipe_K_D
#* _calculate_stagnation_factor
#* _calculate_mean_dw_mass_per_segment
#* _calculate_peak_dw_mass_per_segment
