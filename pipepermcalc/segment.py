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
    ----------
    #ah_todo add attributes

    partitioning_a_dh: float
        Coefficient for correcting the partitioning coefficient for temperature. 
        From regression analysis, a is the slope, see table 5-6 in 
        KWR 2016.056. Constant equal to 7.92169801506708. 
    partitioning_b_dh: float,
        Coefficient for correcting the partitioning coefficient for temperature. 
        From regression analysis, b is the intercept, see table 5-6 in 
        KWR 2016.056. Constant equal to -17.1875608983359. 
    diffusion_a_dh: float 
        Coefficient for correcting the diffusion coefficient for temperature. 
        From regression analysis, a is the slope, see table 5-6 in 
        KWR 2016.056. Constant equal to 61.8565740136974. 
    diffusion_b_dh: float
        Coefficient for correcting the diffusion coefficient for temperature. 
        From regression analysis, b is the intercept, see table 5-6 in 
        KWR 2016.056. Constant equal to -78.9191401984509. 
    assessment_factor_groundwater: integer 
        Factor used to correct calculations for observations in actual pipe 
        permeation. Permeation of PE house connections in groundwater = 3, 
        other pipe materials = 1. See section 7.2 in KWR 2016.056
    assessment_factor_soil: integer
        Factor used to correct calculations for observations in actual pipe 
        permeation. All pipe materials = 1.
    partitioning_a_c: float
        Constant used in the correction for the partitioning coefficent due to 
        the influence of temperature. See equation 5-20 in KWR 2016.056, for 
        partitioning a_c = 0.103965019849463.
    partitioning_Cref_Sw: float
        Reference concentration used in the correction for the partitioning 
        coefficent due to the influence of temperature. Ssee section 5.4.7 in 
        KWR 2016.056. For partitioning, Cref_SW = 1.0.
    diffusion_a_c: float
        Constant used in the correction for the diffusion coefficent due to 
        the influence of temperature. See equation 5-18 in KWR 2016.056, for 
        diffusion a_c = 0.784077209735583.
    diffusion_Cref_Sw: float
        Reference concentration used in the correction for the diffusion 
        coefficent due to the influence of temperature. Ssee section 5.4.6 in 
        KWR 2016.056. For partitioning, Cref_SW = 0.5.

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
        permeation_direction: string #ah_todo enum?? @Bram -> limit choice of direction
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
    volume: float
        Volume of the pipe segment

    reference_pipe_material_dict: dictionary 
        Reference dictionary with regression values (slope=a-value, intercept=b-value)
        for the partitioning (K) and diffusion (D) coefficient for the available 
        plastics (PE40, PE80, SBR, EPDM) per chemical groups. 
        Used to calculate the specifice partitioning or diffusion coefficient 
        from the reference value. 
        Chemical group numbers: Expert opinion (M. Meerkerk), 
        see table 5-4 KWR 2016.056, Group 1: PAK, MAK, ClArom, ClAlk, Arom, Alk
        Group 2: PCB, Group 3: overig, onbekend, O2, Cl, BDE

    log_Kpw_ref: float
        partitioning coefiicient under lab conditions, [-]
    f_Ktemp: float
        Temperature correction factor for partitioning coefficient, [-]
    f_Kconc: float
        Concentration correction factor for partitioning coefficient, [-]
    log_Kpw: float
        Calculated log partitioning coefficient for the given chemical and pipe material, [-]
    log_Dp_ref: float
        Diffusion coefficient under lab conditions, m2/s
    f_Dtemp:float
        Temperature correction factor for diffusion coefficient, [-]
    f_Dconc: float
        Concentration correction factor for diffusion coefficient, [-]
    log_Dp: float
        Calculated log diffusion coefficient for the given chemical and pipe material, [-]
    permeation_coefficient: float
        Calculated permeation coefficient for plastic-water, m2/day
    stagnation_factor: float
        Correction for the decrease in the concentratino gradient near the 
        inner wall of the pipe during stagnation (e.g. no flow at at night)


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
        Creates a pipe segment with the attributes of the pipe (length, 
        thickness, diameter, material etc.). 
        
        Parameters
        ----------
        name: string
            name of the pipe segment
        material: string enum?? @Bram -> set choice of materials
            e.g. PE40, PE80, PVC, EPDM, rubber etc.
        length: float
            Length of pipe segment, meters 
        inner_diameter: float
            Inner diameter of pipe segment, meters
        thickness: float
            Thickness of pipe segment, meters
        permeation_direction: string #ah_todo enum?? @Bram -> limit choice of direction
            Direction of permeation through the pipe segment. Options are 
            'perpendicular' or 'parallel'. Default permeation is perpendicular 
            to the flow direction. See schematic XX in read the docs.
        diffusion_path_length: float
            In the case of permeation perpendicular to the flow direction, a 
            diffusion path length is required to calculate the permeation 
            through the pipe segment. For example in the case of a pipe 
            coupling rings. If no value is given, diffusion is assumed 
            perpendicular to the flow direction and the thickness is 
            used to calculate the diffusion through the pipe segment. 
            Unit meters.

        '''  


        #Constants for various LogK and Log D equations
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

        self.length = float(length)
        self.inner_diameter = float(inner_diameter)
        self.thickness = float(thickness)
        self.permeation_direction = permeation_direction
    
        if diffusion_path_length is None:
            self.diffusion_path_length = self.thickness
        else:
            self.diffusion_path_length = float(diffusion_path_length)    
        
        # Checks here that length, inner diameter, thickness, diffusion path length > 0
        check_values = ['length', 'inner_diameter', 'thickness', 'diffusion_path_length']
        for check_value in check_values:
            value = getattr(self, check_value)
            if value <= 0:
                raise ValueError(f'{check_value} must be > 0 ')

        outer_diameter = inner_diameter + thickness

        if permeation_direction == 'parallel':
            volume = 0 
            permeation_surface_area = (math.pi * ((inner_diameter + thickness) ** 2 - inner_diameter ** 2) )/4
        elif permeation_direction == 'perpendicular':
            volume = math.pi * (inner_diameter / 2) ** 2 * length
            permeation_surface_area =(math.pi * inner_diameter * length)

                   
        self.volume = volume
        self.permeation_surface_area = permeation_surface_area
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter


    # @ah_todo revert back to csv? seperate file? 
    # From Bram, @MartinvdS-> suggest to implement the "named tuple" method, leave for now do at the end
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

        See table 5-3, equations 5-17 and 5-19 in KWR 2016.056

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

        Cg_Sw = min(pipe_permeability_dict['concentration_groundwater'] / pipe_permeability_dict['solubility'], 1)
        f_conc = a_c * (Cg_Sw - Cref_Sw)

        return f_conc


    def _correct_for_age(self,):
        '''
        Age correction, none implemented yet'''

        f_age = 0.000
        return f_age


    def _calculate_ref_logK(self,
                           pipe_permeability_dict):
        '''Calculate the reference log K'''

        a_ref = self.reference_pipe_material_dict[self.material]['ref_log_K_a'][pipe_permeability_dict['chemical_group_number']]
        b_ref = self.reference_pipe_material_dict[self.material]['ref_log_K_b'][pipe_permeability_dict['chemical_group_number']]
        log_Kpw_ref = a_ref * pipe_permeability_dict['log_octanol_water_partitioning_coefficient'] + b_ref

        return log_Kpw_ref


    def _calculate_ref_logD(self,
                           pipe_permeability_dict):
        '''Calculate the reference log K'''

        a_ref = self.reference_pipe_material_dict[self.material]['ref_log_D_a'][pipe_permeability_dict['chemical_group_number']]
        b_ref = self.reference_pipe_material_dict[self.material]['ref_log_D_b'][pipe_permeability_dict['chemical_group_number']]
        log_Dp_ref = a_ref * pipe_permeability_dict['molecular_weight'] + b_ref

        return log_Dp_ref    


    def _calculate_logK(self, 
                        pipe_permeability_dict):
        ''' 
        Calculate the LogK value for the pipe material, correct for temperature,
        concentration and age. Assign the values to the pipe_permeability_dict
        
        See table 5-3 in KWR 2016.056 for explanation of calculations
        
        '''

        # calculate reference log K plastic-water (log kpw) 
        log_Kpw_ref = self._calculate_ref_logK(pipe_permeability_dict=pipe_permeability_dict, )

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
        
        self.log_Kpw_ref = log_Kpw_ref
        self.f_Ktemp = f_Ktemp   
        self.f_Kconc = f_Kconc
        self.log_Kpw = log_Kpw

        return log_Kpw


    def _calculate_logD(self, 
                        pipe_permeability_dict):
        ''' 
        Calculate the LogK value for the pipe material, correct for temperature,
        concentration and age. Assign the values to the pipe_permeability_dict

        See table 5-3 in KWR 2016.056 for explanation of calculations
        
        '''
        
        # calculate reference log D plastic (log Dp) 
        log_Dp_ref = self._calculate_ref_logD(pipe_permeability_dict=pipe_permeability_dict, )

        # correct for temperature, concentration, age
        f_Dtemp = self._correct_for_temperature(pipe_permeability_dict=pipe_permeability_dict, 
                        coefficient_name ='molecular_weight', 
                            a_dh = self._diffusion_a_dh, 
                            b_dh = self._diffusion_b_dh, 
                        )

        f_Dconc = self._concentration_correction(pipe_permeability_dict=pipe_permeability_dict,
                                a_c = self.diffusion_a_c , 
                                Cref_Sw = self.diffusion_Cref_Sw) 
        
        f_Dage = self._correct_for_age()

        # sum corrections for final Log D
        log_Dp = log_Dp_ref + f_Dtemp + f_Dconc + f_Dage



        #ah discussed wtih Bram not to store these values, only calculate them on the fly
        self.log_Dp_ref = log_Dp_ref
        self.f_Dtemp = f_Dtemp    
        self.f_Dconc = f_Dconc
        self.log_Dp = log_Dp #m2/s

        return log_Dp

    def _calculate_permeation_coefficient(self,):
        ''' Calculate the permeation coefficient for the segment'''
        #Permeation coefficient for plastic-water (Ppw), unit: m2/day
        permeation_coefficient = (24 * 60 * 60 * 
                                        (10 ** self.log_Dp) 
                                        * 10 ** self.log_Kpw)
        return permeation_coefficient

    def _calculate_pipe_K_D(self,
                            pipe_permeability_dict,
                            _groundwater_conditions_set, 
        ):
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
            self.log_Kpw = self._calculate_logK(pipe_permeability_dict = pipe_permeability_dict)

            # calculate log D plastic (log Dp) 
            self.log_Dp = self._calculate_logD(pipe_permeability_dict = pipe_permeability_dict)

            #Permeation coefficient for plastic-water (Ppw), unit: m2/day
            self.permeation_coefficient = self._calculate_permeation_coefficient()
            

    def _calculate_stagnation_factor(self,):
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

        stagnation_factor = 10 ** max((((self.log_Dp + 12.5) / 2 + 
                                self.log_Kpw) * 0.73611 + 
                                -1.03574 ), 0)
        return stagnation_factor
    

    def _calculate_mean_dw_mass_per_segment(self, 
                                            pipe_permeability_dict,
                                            concentration_drinking_water,
                                            _groundwater_conditions_set,
                                            flow_rate=None,
                                        ): 
        '''
        Calculates the mean mass in drinking water for a 24 hour period given a 
        groundwater concentration, for each pipe segment.
        Mean mass in drinking water added to the pipe_permeability_dict.
        
        Parameters
        ----------
        '''

        concentration_groundwater = pipe_permeability_dict['concentration_groundwater'] 
         
        # From equation 4-7 in KWR 2016.056, but not simplifying the mass flux 
        # in equation 4-5 
        delta_c = concentration_groundwater - concentration_drinking_water

        self.mass_chemical_drinkwater = ((self.permeation_coefficient 
                                          * self.permeation_surface_area 
                                          * delta_c / self.diffusion_path_length ) 
                                            / self.assessment_factor_groundwater)


    def _calculate_peak_dw_mass_per_segment(self, 
                                        pipe_permeability_dict, 
                                        concentration_drinking_water,
                                        _groundwater_conditions_set,
                                        stagnation_time_hours = 8, 
                                        flow_rate=None,
                                        ):
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

        concentration_groundwater = pipe_permeability_dict['concentration_groundwater'] 

        self.stagnation_factor = self._calculate_stagnation_factor()
        delta_c = concentration_groundwater - concentration_drinking_water

        # From equation 4-10 KWR 2016.056, but not simplifying the mass flux 
        # in equation 4-5 and rearranging to remove C_dw from the equation       
        self.mass_chemical_drinkwater = ((self.permeation_coefficient 
                                             * self.permeation_surface_area 
                                             * delta_c / self.diffusion_path_length 
                                             * stagnation_time * self.stagnation_factor) 
                                            / self.assessment_factor_groundwater)


    def _update_partitioning_coefficient(self, 
                                        new_log_Kpw=None,):
        ''' Function to update the partitioning coefficient and the associated 
        permeation coefficient

        Parameters
        ----------
        new_log_Kpw:float
            New value for the partitioning coefficient for the given pipe segment
        segment_name: string
            name of the pipe segment        
        '''

        self.log_Kpw = new_log_Kpw

        #Permeation coefficient for plastic-water (Ppw), unit: m2/day
        self.permeation_coefficient = self._calculate_permeation_coefficient()
       
        
    def _update_diffusion_coefficient(self, 
                                        new_log_Dp=None, ):
        ''' Function to update the diffusion coefficient and the associated 
        permeation coefficient

        Parameters
        ----------
        new_log_Dp:float
            New value for the diffusion coefficient for the given pipe segment
        segment_name: string
            name of the pipe segment        
        '''

        self.log_Dp = new_log_Dp

        #Permeation coefficient for plastic-water (Ppw), unit: m2/day
        self.permeation_coefficient = self._calculate_permeation_coefficient()

