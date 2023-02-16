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

class Segment:
    ''' 
    Segment object class to make segments of the pipe.

    Attributes
    -------
    #ah_todo add attributes
    name: string
        name of the pipe segment
    material: enum?? @Bram -> set choice of materials
        e.g. PE40, PE80, PVC, EPDM, rubber etc.
    length: float
        Length of pipe segment, meters 
    inner_diameter: float
        Inner diameter of pipe segment, meters
    thickness: float
        Thickness of pipe segment, meters
    permeation_direction: string
        Direction of permeation through the pipe segment. Options are 
        'perpendicular' or 'parallel'. Default permeation is perpendicular 
        to the flow direction. See schematic XX in read the docs.
        @AH_todo add schematic of permeation for parallel vs. perpendicular
    diffusion_path_length: float
        In the case of permeation perpendicular to the flow direction, a 
        diffusion path length is required to calculate the permeation 
        through the pipe segment. For example in the case of a pipe 
        coupling rings. If no value is given, diffusion is assumed 
        perpendicular to the flow direction and the thickness is 
        used to calculate the diffusion through the pipe segment. 
        Unit meters.
    '''

    count = 0 # count of pipe segments

    def __init__(self, 
                name=None,
                material=None,
                length=None,
                inner_diameter=None,
                thickness=None,
                permeation_direction='perpendicular',
                diffusion_path_length=None,
                ):
        '''
        Creates a pipe segment. Creates pipe_dictionary as attribute 
        containing information on the different pipe segment materials and sizes. 
        #ah_todo update these definitions

        '''  
        #constants
        self._partitioning_a_dh = 7.92169801506708 #see table 5-6 in KWR 2016.056
        self._partitioning_b_dh = -17.1875608983359 #see table 5-6 in KWR 2016.056
        self._diffusion_a_dh = 61.8565740136974 #see table 5-6 in KWR 2016.056
        self._diffusion_b_dh = -78.9191401984509 #see table 5-6 in KWR 2016.056
        self.assessment_factor_groundwater = 3 
        self.assessment_factor_soil = 1
        self.partitioning_a_c = 0.103965019849463 #see equation 5-20 in KWR 2016.056
        self.partitioning_Cref_Sw = 1.000 #see section 5.4.7 in KWR 2016.056
        self.diffusion_a_c = 0.784077209735583 #see equation 5-18 in KWR 2016.056
        self.diffusion_Cref_Sw = 0.5 #see section 5.4.6 in KWR 2016.056

        self.name = name
        self.material = material
        self.length = length
        self.inner_diameter = inner_diameter
        self.thickness = thickness
        self.permeation_direction = permeation_direction
        self.diffusion_path_length = diffusion_path_length    

        if diffusion_path_length is None:
            diffusion_path_length = thickness
        else:
            pass
        
        outer_diameter = inner_diameter + thickness

        if permeation_direction == 'parallel':
            volume = 0 
            # outer_surface_area = None #not applicable to this type of permeation
            permeation_surface_area = (math.pi * ((inner_diameter + thickness) ** 2 - inner_diameter ** 2) )/4
        elif permeation_direction == 'perpendicular':
            volume = math.pi * (inner_diameter / 2) ** 2 * length
            # outer_surface_area = (math.pi * outer_diameter * length)
            # inner_surface_area = (math.pi * inner_diameter * length)
            permeation_surface_area =(math.pi * inner_diameter * length)

        # @MartinvdS check about drawing #1 = outer diameter used, 
        # #3 inner diameter used for SA calculations, see notebook and 
        # sheet "dimensies tertiare"

        pipe_dictionary = {
                # 'number_segments': self.count,
                # 'segment_list': self.segment_list,
                'total_length':length,
                'total_volume':volume,
                'material': material,
                'length': length,
                'outer_diameter': outer_diameter,
                'inner_diameter': inner_diameter,
                'thickness': thickness,
                'diffusion_path_length': diffusion_path_length,
                'volume': volume,
                # 'outer_surface_area': outer_surface_area,
                # 'inner_surface_area': inner_surface_area,
                'permeation_surface_area': permeation_surface_area,
                'permeation_direction':permeation_direction,
                }
            
        self.pipe_dictionary = pipe_dictionary


    # @ah_todo revert back to csv? seperate file? 
    # @MartinvdS-> suggest to implement the "named tuple" method, leave for now do at the end
    reference_pipe_material_dict = \
        {
        "PE40": {
            "ref_log_D_a": {
                1: -0.011,	
                2: -0.00629,
                3: -0.006,
                },
            "ref_log_D_b": {
                1: -10.688,
                2: -11.000,
                3: -11.000
                },
            "ref_log_K_a": {
                1: 1.097,
                2: 1.059,
                3: 0.979
                },
            "ref_log_K_b": {
                1: -0.689,
                2: -0.67,
                3: -1.796,
                },
        },
        "PE80": {
            "ref_log_D_a": {
                1: -0.011,	
                2: -0.00629,
                3: -0.00629,
                },
            "ref_log_D_b": {
                1: -11.188,
                2: -11.188,
                3: -11.500,
                },
            "ref_log_K_a": {
                1: 1.185,	
                2: 1.185,
                3: 1.231,
                },
            "ref_log_K_b": {
                1: -1.437,	
                2: -1.437,
                3: -2.606,
                },
        },
        "SBR": {
            "ref_log_D_a": {
                1: 0.950647410867427,	
                2: 0.950647410867427,
                3: 0.950647410867427,
                },
            "ref_log_D_b": {
                1: 0.0,
                2: 0.0,
                3: 0.0,
                },
            "ref_log_K_a": {
                1: 1.0452,	
                2: 1.0452,
                3: 1.0452,
                },
            "ref_log_K_b": {
                1: -0.3686,	
                2: -0.3686,
                3: -0.3686,
                },
        },  
        "EPDM": {
            "ref_log_D_a": {
                1: 0.920996123470591,	
                2: 0.920996123470591,
                3: 0.920996123470591,
                },
            "ref_log_D_b": {
                1: 0.0,
                2: 0.0,
                3: 0.0,
                },
            "ref_log_K_a": {
                1: 1.0675,	
                2: 1.0675,
                3: 1.0675,
                },
            "ref_log_K_b": {
                1: -0.3002,	
                2: -0.3002,
                3: -0.3002,
                },
        },      
        }

    def _correct_for_temperature(self,
                                pipe_permeability_dict=None, 
                                coefficient_name=None, 
                                a_dh=None,
                                b_dh=None,
                                ):
        '''
        Temperature correction for the partitioning and diffusion coefficients, 
        
        See table 5-3 in KWR 2016.056

        Parameters
        ----------
        pipe_permeability_dict: dictionary
            Dictionary of permeability coefficients
        temperature_groundwater: float
            Temperature of the groundwater, degrees Celcius
        coefficient_name: string
            Either 'solubility' or 'molecular_weight'
        a_dh: string
            Coefficient for correcting the partitioning or diffusion coefficient
        b_dh: string
            Coefficient for correcting the partitioning or diffusion coefficient

        Returns
        -------
        f_temp: string
            Temperature correction factor for the partitioning or diffusion 
            coefficient
        '''
        temperature_groundwater=pipe_permeability_dict['temperature_groundwater']

        R = 0.008314 #universal gas constant
        reference_temperature = 25 # deg. C
        dh = a_dh * math.log10(pipe_permeability_dict[coefficient_name]) + b_dh
        f_temp = dh / (R * math.log(10)) * (1 / (reference_temperature + 273) - 1 / (temperature_groundwater + 273))
        return f_temp


    def _concentration_correction(self,
                         pipe_permeability_dict=None, 
                        a_c=None,
                        Cref_Sw=None):
        '''
        Correction factor for the influence of concentration on the 
        partitioning or diffusion coefficient 

        See table 5-3 in KWR 2016.056

        Parameters
        ----------
        pipe_permeability_dict: dictionary
            Dictionary of permeability coefficients

        Returns
        -------
        f_conc: string
            Concentration correction factor for the partitioning or diffusion 
            coefficient
        '''

        # from ppc_database material column K27-29
        Cg_Sw = min(pipe_permeability_dict['concentration_groundwater'] / pipe_permeability_dict['solubility'], 1)
        f_conc = a_c * (Cg_Sw - Cref_Sw)
        return f_conc


    def _correct_for_age(self,):
        '''
        Age correction, none implemented yet'''

        f_age = 0.000
        return f_age
    

    def _calculate_logK(self, 
                        pipe_permeability_dict):
        ''' 
        Calculate the LogK value for the pipe material, correct for temperature,
        concentration and age. Assign the values to the pipe_permeability_dict
        
        See table 5-3 in KWR 2016.056 for explanation of calculations
        
        Parameters
        ----------
        pipe_material: string
            Choice of pipe material: PE40, PE80, SBR, EPDM
        '''

        # calculate reference log K plastic-water (log kpw) 
        a_ref = self.reference_pipe_material_dict[self.material]['ref_log_K_a'][pipe_permeability_dict['chemical_group_number']]
        b_ref = self.reference_pipe_material_dict[self.material]['ref_log_K_b'][pipe_permeability_dict['chemical_group_number']]
        log_Kpw_ref = a_ref * pipe_permeability_dict['log_octanol_water_partitioning_coefficient'] + b_ref

        # correct for temperature, concentration, age
        f_Ktemp = self._correct_for_temperature(pipe_permeability_dict=pipe_permeability_dict, 
                        coefficient_name = 'solubility',
                            a_dh = self._partitioning_a_dh, 
                            b_dh = self._partitioning_b_dh, 
                        )

        f_Kconc = self._concentration_correction(pipe_permeability_dict=pipe_permeability_dict,
                                a_c = self.partitioning_a_c,
                                Cref_Sw = self.partitioning_Cref_Sw) 
        
        f_Kage = self._correct_for_age()

        # sum corrections for final Log k
        log_Kpw = log_Kpw_ref + f_Ktemp + f_Kconc + f_Kage

        #ah discussed wtih Bram not to store these values, only calculate them on the fly
        # pipe_permeability_dict['log_Kpw_ref'] = log_Kpw_ref
        # pipe_permeability_dict['f_Ktemp'] = f_Ktemp   
        # pipe_permeability_dict['f_Kconc'] = f_Kconc
        # pipe_permeability_dict['log_Kpw'] = log_Kpw

        return log_Kpw


    def _calculate_logD(self, 
                        pipe_permeability_dict):
        ''' 
        Calculate the LogK value for the pipe material, correct for temperature,
        concentration and age. Assign the values to the pipe_permeability_dict

        See table 5-3 in KWR 2016.056 for explanation of calculations
        
        Parameters
        ----------
        pipe_material: string
            Choice of pipe material: PE40, PE80, SBR, EPDM
        '''
        
        # calculate reference log D plastic (log Dp) 
        a_ref = self.reference_pipe_material_dict[self.material]['ref_log_D_a'][self.pipe_permeability_dict['chemical_group_number']]
        b_ref = self.reference_pipe_material_dict[self.material]['ref_log_D_b'][self.pipe_permeability_dict['chemical_group_number']]
        log_Dp_ref = a_ref * self.pipe_permeability_dict['molecular_weight'] + b_ref

        # correct for temperature, concentration, age
        f_Dtemp = self._correct_for_temperature(pipe_permeability_dict=self.pipe_permeability_dict, 
                        temperature_groundwater=self.temperature_groundwater, 
                        coefficient_name ='molecular_weight', 
                            a_dh = self._diffusion_a_dh, 
                            b_dh = self._diffusion_b_dh, 
                        )

        f_Dconc = self._concentration_correction(pipe_permeability_dict=self.pipe_permeability_dict,
                                a_c = self.diffusion_a_c , 
                                Cref_Sw = self.diffusion_Cref_Sw) 
        
        f_Dage = self._correct_for_age()

        # sum corrections for final Log D
        log_Dp = log_Dp_ref + f_Dtemp + f_Dconc + f_Dage

        segment_dict['log_Dp_ref'] = log_Dp_ref
        segment_dict['f_Dtemp'] = f_Dtemp    
        segment_dict['f_Dconc'] = f_Dconc
        segment_dict['log_Dp'] = log_Dp #m2/s

        return segment_dict


    def _calculate_pipe_K_D(self,
                            pipe_permeability_dict,
                            _groundwater_conditions_set, 
                                 chemical_name, 
                                 temperature_groundwater, 
                                 concentration_groundwater): #ah_funct1
        '''
        Fetch the pipe and chemical information corresponding to the given pipe 
        material and chemical choice. Creates the pipe_permeability_dict, 
        which consists of chemical and permeability related coefficients for the
        different pipe segments.

        See table 5-3 in KWR 2016.056 for explanation of calculations

        Parameters
        ----------
        pipe_material: string
            Choice of pipe material: PE40, PE80, SBR, EPDM
        temperature_groundwater: float
            Temperature of groundwater, in degrees Celcius
        concentration_groundwater: float
            concentration of given chemical in groundwater, in mg/L
        '''

        # Check if the groundwater conditions have been set, if not raise error
        if _groundwater_conditions_set is None:
            raise ValueError('Error, groundwater conditions have not been set. \
                             To set groundwater conditions use .set_groundwater_conditions()')
        else:           

            # calculate log K plastic-water (log kpw) 
            self._calculate_logK(pipe_permeability_dict = pipe_permeability_dict)

            # calculate log D plastic (log Dp) 
            self._calculate_logD(pipe_permeability_dict = pipe_permeability_dict)

            #Permeation coefficient for plastic-water (Ppw), unit: m2/day
            segment_dict['permeation_coefficient'] = (24 * 60 * 60 * 
                                        (10 ** segment_dict['log_Dp']) 
                                        * 10 ** segment_dict['log_Kpw'])
            
            self.pipe_permeability_dict['segments'][segment_name] = segment_dict        

    def _calculate_stagnation_factor(self, 
                                     pipe_segment=None):
        ''' Calculates the stagnation factor given a pipe segment
        
        Parameters
        ----------
        pipe_segment: string
            name of the pipe segment        
        
        Returns
        ----------
        stagnation_factor: string
            Correction for the decrease in the concentratino gradient near the 
            inner wall of the pipe during stagnation (e.g. no flow at at night)

        '''

        stagnation_factor = 10 ** max((((self.pipe_permeability_dict['segments'][pipe_segment]['log_Dp'] + 12.5) / 2 + 
                                self.pipe_permeability_dict['segments'][pipe_segment]['log_Kpw']) * 0.73611 + 
                                -1.03574 ), 0)
        return stagnation_factor
    

    def _calculate_mean_dw_mass_per_segment(self, 
                                            pipe_permeability_dict,
                                            _groundwater_conditions_set,
                                            flow_rate=None,
                                            assessment_factor_groundwater=None,
                                        ): #ah_funct1
        '''
        Calculates the mean mass in drinking water for a 24 hour period given a 
        groundwater concentration, for each pipe segment.
        Mean mass in drinking water added to the pipe_permeability_dict.
        
        Parameters
        ----------
        '''
        segment_surface_area = self.pipe_dictionary['permeation_surface_area']
        segment_volume = self.pipe_dictionary['volume']
        segment_diffusion_path_length = self.pipe_dictionary['diffusion_path_length']         
        segment_diffusion_path_length = self.pipe_dictionary['diffusion_path_length']
        concentration_groundwater = pipe_permeability_dict['concentration_groundwater'] 

        self._calculate_pipe_K_D(pipe_permeability_dict, 
                                 _groundwater_conditions_set, 
                            ) #ah_funct1

        permeation_coefficient = self.pipe_dictionary['permeation_coefficient']

        #@martinvdS, but if we assign some volumes to zero, this messes this up, 
        # so input the actual volume calculation instead?     
           
        # From equation 4-7 in KWR 2016.056, but not simplifying the mass flux 
        # in equation 4-5 and rearranging to remove C_dw from the equation
        mass_drinkwater = ((concentration_groundwater * segment_volume * 
                           segment_surface_area * permeation_coefficient) / 
                            (segment_diffusion_path_length * assessment_factor_groundwater *
                              flow_rate + permeation_coefficient * segment_surface_area))
        
        self.self.pipe_dictionary['mass_drinkwater'] = mass_drinkwater #ah_todo rename to chemical_mass_drinkingwater

    def _calculate_peak_dw_mass_per_segment(self, 
                                         pipe_segment=None,
                                        stagnation_time_hours = 8, ):
        #Segment class() ah_todo: move all functions on segments to the segment class
        '''
        Calculates the peak (maximum) mass in drinking water for a 
        given a stagnation period given a groundwater concentration, for each pipe segment.
        Stagnation period default of 8 hours. Peak mass in drinking 
        water added to the pipe_permeability_dict.
        
        Parameters
        ----------
        pipe_segment: string
            name of the pipe segment        
        stagnation_time_hours: float
            time in hours, default 8 hours

        '''
        stagnation_time = stagnation_time_hours / 24 # days
        segment_volume = self.pipe_dictionary['segments'][pipe_segment]['volume']
        segment_surface_area = self.pipe_dictionary['segments'][pipe_segment]['permeation_surface_area']

        segment_diffusion_path_length = self.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 
        concentration_groundwater = self.pipe_permeability_dict['concentration_groundwater']
        segment_diffusion_path_length = self.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length']
        permeation_coefficient = self.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient']
        stagnation_factor = self._calculate_stagnation_factor(pipe_segment=pipe_segment)

        # From equation 4-10 KWR 2016.056, but not simplifying the mass flux 
        # in equation 4-5 and rearranging to remove C_dw from the equation
        mass_drinkwater =(( permeation_coefficient * segment_surface_area * 
                                stagnation_time * concentration_groundwater *  segment_volume) / 
                                ((segment_diffusion_path_length * self.assessment_factor_groundwater * 
                                  stagnation_factor * segment_volume) + 
                                  (permeation_coefficient * segment_surface_area * stagnation_time) )  )
       
        self.pipe_permeability_dict['segments'][pipe_segment]['mass_drinkwater'] = mass_drinkwater
        self.pipe_permeability_dict['segments'][pipe_segment]['stagnation_factor'] = stagnation_factor

