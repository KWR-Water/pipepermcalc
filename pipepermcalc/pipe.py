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

class Pipe:
    '''
    Pipe object class to make segments of the pipe and calculate the peak and
    mean concentration of a chemical in groundwater and soil.

    Attributes
    -------
    segment_list: LIST,
        list of the names of the different pipe sements

    count: integer
        Count of the number of segments created for the pipe

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
    diffusion_b_dh: float,
        Coefficient for correcting the diffusion coefficient for temperature. 
        From regression analysis, b is the intercept, see table 5-6 in 
        KWR 2016.056. Constant equal to -78.9191401984509. 
    assessment_factor_groundwater: integer, 
        Factor used to correct calculations for observations in actual pipe 
        permeation. Permeation of PE house connections in groundwater = 3, 
        other pipe materials = 1. See section 7.2 in KWR 2016.056
    assessment_factor_soil: integer,
        Factor used to correct calculations for observations in actual pipe 
        permeation. All pipe materials = 1.
    partitioning_a_c: float,
        Constant used in the correction for the partitioning coefficent due to 
        the influence of temperature. See equation 5-20 in KWR 2016.056, for 
        partitioning a_c = 0.103965019849463.
    partitioning_Cref_Sw: float,
        Reference concentration used in the correction for the partitioning 
        coefficent due to the influence of temperature. Ssee section 5.4.7 in 
        KWR 2016.056. For partitioning, Cref_SW = 1.0.
    diffusion_a_c: float,
        Constant used in the correction for the diffusion coefficent due to 
        the influence of temperature. See equation 5-18 in KWR 2016.056, for 
        diffusion a_c = 0.784077209735583.
    diffusion_Cref_Sw: float,
        Reference concentration used in the correction for the diffusion 
        coefficent due to the influence of temperature. Ssee section 5.4.6 in 
        KWR 2016.056. For partitioning, Cref_SW = 0.5.

    reference_pipe_material_dict: dictionary 
        Reference dictionary with regression values (slope=a-value, intercept=b-value)
        for the partitioning (K) and diffusion (D) coefficient for the available 
        plastics (PE40, PE80, SBR, EPDM) per chemical groups. 
        Used to calculate the specifice partitioning or diffusion coefficient 
        from the reference value. 
        Chemical group numbers: Expert opinion (Martin Meerkerk), 
        see table 5-4 KWR 2016.056, Group 1: PAK, MAK, ClArom, ClAlk, Arom, Alk
        Group 2: PCB, Group 3: overig, onbekend, O2, Cl, BDE

    pipe_dictionary: dictionary
        number_segments: float,
            total number of pipe segments
        segment_list: list,
            list of the names of the different pipe sements
        total_length: float,
            sum of the lengths of the pipe segments
        total_volume: float,
            sum of the volumes of the pipe segments
        total_outer_surface_area: float,
            sum of the surface area of the pipe segments
        flow_rate: float,
            flow rate in the pipe, m3/day 
        segments: dictionary
            dictionary of the individual pipe segments, containing the 
            segment material, segment length (m), outer outer_diameter (m), thickness (m), 
            volume (m3) and surface area (m2)
    
    pipe_permeability_dict: dictionary
        CAS_number: string,
            CAS is a unique identification number assigned by the Chemical 
            Abstracts Service (CAS)
        chemical_name_NL: string,
            Name of the chemical given in Dutch
        chemical_name: string,
            @Bram to be replaced by the function restricting input of names
        molecular_weight: float,
            Mass of one mole of a given chemical, g/mol
        solubility:	float,
            solubility of given chemical in water, g/m3
        log_octanol_water_partitioning_coefficient:	float,
            Partition coefficient for the two-phase system consisting of 
            n-octanol and water, Log Kow, [-]
        log_distribution_coefficient: float,
            Ratio of the amount of chemical  adsorbed onto soil per amount 
            of water, m3/g
        chemical_group: string,
            Grouping of chemicals (expert opinion) with similar properties
            for permeation: Group 1: PAK, MAK, ClArom, ClAlk, Arom, Alk, 
            Group 2: PCB, Group 3: overig, onbekend, O2, Cl, BDE. See KWR 2016.056
        chemical_group_number: integer,
            Integer corresponding to the chemical group 
        molecular_volume: float,
            Volume occupied by one mole of the substance at a given 
            temperature and pressure, cm3/mol
        Drinking_water_norm: float,
            Concentration allowable in the Dutch Drinking water decree, ug/L
        concentration_groundwater: float,
            Concentration of the given chemical in groundwater, g/m3
        temperature_groundwater: float,
            Temperature of the groundwater, degrees Celcius
        log_Kpw_ref: float,
            partitioning coefiicient under lab conditions, [-]
        f_Ktemp: float,
            Temperature correction factor for partitioning coefficient, [-]
        f_Kconc: float,
            Concentration correction factor for partitioning coefficient, [-]
        log_Kpw: float,
            Calculated log partitioning coefficient for the given chemical and pipe material, [-]
        log_Dp_ref: float,
            Diffusion coefficient under lab conditions, m2/s
        f_Dtemp:float,
            Temperature correction factor for diffusion coefficient, [-]
        f_Dconc: float,
            Concentration correction factor for diffusion coefficient, [-]
        log_Dp: float,
            Calculated log diffusion coefficient for the given chemical and pipe material, [-]
        permeation_coefficient: float,
            Calculated permeation coefficient for plastic-water, m2/day
        stagnation_time_hours:float,
            Time in hours which water in pipe is stagnant, hours.
        flux_max_stagnation: float,
            Maximum flux during stagnation period to remain below the drinking 
            water standard, mg/day
        flux_max_stagnation_per_m2: float,
            Maximum flux per square meter (surface area) of pipe during 
            stagnation period to remain below the drinking water standard, g/day
        flux_max_per_day: float,
            Maximum flux in a day (24 hours) to remain below the drinking 
            water standard, g/day
        flux_max_per_day_per_m2: float,
            Maximum flux per square meter (surface area) of pipe in one day 
            (24 hours) to remain below the drinking water standard, g/day        
        stagnation_factor: float,
            Correction for the decrease in the concentratino gradient near the 
            inner wall of the pipe during stagnation (e.g. no flow at at night)
        concentration_peak_without_stagnation: float,
            Concentration in groundwater which, wihtout a stagnation period, 
            would not result in a peak concentration in drinking water exceeding 
            the drinking water norm, g/m3
        concentration_peak_after_stagnation: float, 
            Concentration in groundwater which, after a stagnation period, 
            would not result in a peak concentration in drinking water exceeding 
            the drinking water norm, g/m3
        concentration_peak_soil: float,
            Concentration in soil which, after a stagnation period, 
            would not result in a peak concentration in drinking water exceeding 
            the drinking water norm, mg/kg
        concentration_mean: float,
            Mean concentration in groundwater which would would not result in 
            a mean daily (24 horus) concentration in drinking water exceeding 
            the drinking water norm, g/m3
        concentration_mean_soil: float,      
            Mean concentration in soil which would would not result in 
            a mean daily (24 horus) concentration in drinking water exceeding 
            the drinking water norm, mg/kg
    '''

    count = 0 # count of pipe segments

    def __init__(self, 
                ):
        '''
        
        '''       
        segment_list = [] # list of pipe segment names
        self.segment_list = segment_list
        self._groundwater_conditions_set = False
        self._flow_rate_set = False

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

       
    # @ah_todo revert back to csv, seperate file @Bram will think about this
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


    def set_groundwater_conditions(self,
                                   chemical_name=None,                                    
                                   concentration_groundwater=None,
                                   temperature_groundwater=None):
        ''' 
        Specifies the chemical of interest, concentration and temperature in the 
        groundwater and returns the parameters as attributes of the class.
        
        Parameters
        ----------
        chemical_name: string, 
            @Bram to be replaced by the function restricting input of names
        concentration_groundwater: float,
            Concentration of the given chemical in groundwater, g/m3
        temperature_groundwater: float,
            Temperature of the groundwater, degrees Celcius
        '''
        
        self.concentration_groundwater = concentration_groundwater
        self.temperature_groundwater = temperature_groundwater
        self.chemical_name = chemical_name
        self._groundwater_conditions_set = True
    

    def set_flow_rate(self, 
                      flow_rate=0.5): 
        ''' set this in function by itself
        , also add check, same as groundwater conditions
        
        Parameters
        ----------        
        flow_rate: float
            flow_rate through pipe. Defaul of 0.5 m3/day.
        '''
        self.flow_rate = flow_rate
        self._flow_rate_set = True


    def _fetch_chemical_database(self,
                                chemical_name=None,):
        ''' 
        Fetch the pipe and chemical information corresponding to the given 
        pipe material and chemical choice and creates a dictionary 
        pipe_permeability_dict which consists of chemical and permeability 
        related coefficients.

        Parameters
        ----------
        chemical_name: string 
            @Bram to be replaced by the function restricting input of names
        '''
        
        ppc_database = read_csv(file_path.parent / 'database' / 'ppc_database.csv',  skiprows=[1, 2] ) 

        df = ppc_database[ppc_database['chemical_name'].str.contains(chemical_name)]

        chemical_dict = df.to_dict('records')[0]

        return chemical_dict


    def _correct_for_temperature(self,
                                pipe_permeability_dict=None, 
                                temperature_groundwater=None,
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
        temperature_groundwater: float,
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
        # a_c = 0.103965019849463
        # Cref_Sw = 1.000
        Cg_Sw = min(pipe_permeability_dict['concentration_groundwater'] / pipe_permeability_dict['solubility'], 1)
        f_conc = a_c * (Cg_Sw - Cref_Sw)
        return f_conc


    def _correct_for_age(self,):
        '''
        Age correction, none implemented yet'''

        f_age = 0.000
        return f_age
    

    def _calculate_logK(self,
                        pipe_material=None,
                        segment_dict=None,):
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
        a_ref = self.reference_pipe_material_dict[pipe_material]['ref_log_K_a'][self.pipe_permeability_dict['chemical_group_number']]
        b_ref = self.reference_pipe_material_dict[pipe_material]['ref_log_K_b'][self.pipe_permeability_dict['chemical_group_number']]
        log_Kpw_ref = a_ref * self.pipe_permeability_dict['log_octanol_water_partitioning_coefficient'] + b_ref

        # correct for temperature, concentration, age
        f_Ktemp = self._correct_for_temperature(pipe_permeability_dict=self.pipe_permeability_dict, 
                        temperature_groundwater=self.temperature_groundwater, 
                        coefficient_name = 'solubility',
                            a_dh = self._partitioning_a_dh, 
                            b_dh = self._partitioning_b_dh, 
                        )

        f_Kconc = self._concentration_correction(pipe_permeability_dict=self.pipe_permeability_dict,
                                a_c = self.partitioning_a_c,
                                Cref_Sw = self.partitioning_Cref_Sw) 
        
        f_Kage = self._correct_for_age()

        # sum corrections for final Log k
        log_Kpw = log_Kpw_ref + f_Ktemp + f_Kconc + f_Kage

        segment_dict['log_Kpw_ref'] = log_Kpw_ref
        segment_dict['f_Ktemp'] = f_Ktemp   
        segment_dict['f_Kconc'] = f_Kconc
        segment_dict['log_Kpw'] = log_Kpw

        return segment_dict

    def _calculate_logD(self, 
                       pipe_material=None,
                        segment_dict=None,):
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
        a_ref = self.reference_pipe_material_dict[pipe_material]['ref_log_D_a'][self.pipe_permeability_dict['chemical_group_number']]
        b_ref = self.reference_pipe_material_dict[pipe_material]['ref_log_D_b'][self.pipe_permeability_dict['chemical_group_number']]
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
                            pipe_material=None, 
                            segment_name=None): 
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
        if self._groundwater_conditions_set is None:
            raise ValueError('Error, groundwater conditions have not been set. \
                             To set groundwater conditions use .set_groundwater_conditions()')
        else:           

            if self.count <= 1:
                self.pipe_permeability_dict = self._fetch_chemical_database(chemical_name=self.chemical_name)

                self.pipe_permeability_dict['concentration_groundwater'] = self.concentration_groundwater
                self.pipe_permeability_dict['temperature_groundwater'] = self.temperature_groundwater
                pipe_material = self.pipe_dictionary['segments'][segment_name]['material']
                
                self.pipe_permeability_dict['segments'] = {}

            else:
                pass
            
            segment_dict = {}

            # calculate log K plastic-water (log kpw) 
            segment_dict.update(self._calculate_logK(pipe_material=pipe_material, 
                                    segment_dict=segment_dict))

            # calculate log D plastic (log Dp) 
            segment_dict.update(self._calculate_logD(pipe_material=pipe_material,
                                                        segment_dict=segment_dict))

            #Permeation coefficient for plastic-water (Ppw), unit: m2/day
            segment_dict['permeation_coefficient'] = (24 * 60 * 60 * 
                                        (10 ** segment_dict['log_Dp']) 
                                        * 10 ** segment_dict['log_Kpw'])
            
            self.pipe_permeability_dict['segments'][segment_name] = segment_dict

    #@MartinvdS added this function since the permeation_coefficient also needs 
    # to be updated, otherwise the Log Kpw is only ~half used
    def _update_partitioning_coefficient(self, 
                                        new_log_Kpw=None, 
                                        segment_name=None):
        ''' Function to update the partitioning coefficient and the associated 
        permeation coefficient

        Parameters
        ----------
        new_log_Kpw:float
            New value for the partitioning coefficient for the given pipe segment
        segment_name: string
            name of the pipe segment        
        '''

        self.pipe_permeability_dict['segments'][segment_name]['log_Kpw'] = new_log_Kpw

        #Permeation coefficient for plastic-water (Ppw), unit: m2/day
        self.pipe_permeability_dict['segments'][segment_name]['permeation_coefficient'] = (24 * 60 * 60 * 
                                    (10 ** self.pipe_permeability_dict['segments'][segment_name]['log_Dp']) 
                                    * 10 ** self.pipe_permeability_dict['segments'][segment_name]['log_Kpw'])
       
        
    def _update_diffusion_coefficient(self, 
                                        new_log_Dp=None, 
                                        segment_name=None):
        ''' Function to update the diffusion coefficient and the associated 
        permeation coefficient

        Parameters
        ----------
        new_log_Dp:float
            New value for the diffusion coefficient for the given pipe segment
        segment_name: string
            name of the pipe segment        
        '''

        self.pipe_permeability_dict['segments'][segment_name]['log_Dp'] = new_log_Dp

        #Permeation coefficient for plastic-water (Ppw), unit: m2/day
        self.pipe_permeability_dict['segments'][segment_name]['permeation_coefficient'] = (24 * 60 * 60 * 
                                    (10 ** self.pipe_permeability_dict['segments'][segment_name]['log_Dp']) 
                                    * 10 ** self.pipe_permeability_dict['segments'][segment_name]['log_Kpw'])

    def add_segment(self,
                    name=None,
                    material=None,
                    length=None,
                    outer_diameter=None,
                    thickness=None,
                    diffusion_path_length=None,
                    ):
        
        '''
        Adds a segment to the pipe. Creates pipe_dictionary as attribute 
        containing information on the different pipe segment materials and sizes. 

        Parameters
        ----------
        name: string
            name of the pipe segment
        material: enum?? @Bram -> set choice of materials
            e.g. PE40, PE80, PVC, EPDM, rubber etc.
        length: float
            Length of pipe segment, meters 
        outer_diameter: float
            Outer outer_diameter of pipe segment, meters
            
            @MartinvdS is this the inner or the outer outer_diameter? 

            @Martin, why in excel is the thickness of the material (column H, 
            sheet "dimensies...") the thickness of "both" sides of the material (2 * thickness)
            Also why is the inner outer_diameter used to calculate the dienstkraan 
            and koppelstuk, shouldn't it be the outer outer_diameter?

        thickness: float
            Thickness of pipe segment, meters
        diffusion_path_length: float
            In the case of permeation perpendicular to the flow direction, a 
            diffusion path length is required to calculate the permeation 
            through the pipe segment. For example in the case of a pipe 
            coupling rings. If no value is given, diffusion is assumed 
            perpendicular to the flow direction and the thickness is 
            used to calculate the diffusion through the pipe segment. 
            Unit meters.
        '''

        self.count += 1 #count the number of segments created
        self.segment_list.append(name) 

        if diffusion_path_length is None:
            diffusion_path_length = thickness
        else:
            pass

        inner_diameter = outer_diameter - thickness
        volume = math.pi * (inner_diameter / 2) ** 2 * length
        outer_surface_area = (math.pi * outer_diameter * length)
        inner_surface_area = (math.pi * inner_diameter * length)
        # @MartinvdS check about the contact area used, 
        # see sheet 'dimensies tertiare structuur', column K

        if self.count >1:

            total_length = length + self.pipe_dictionary['total_length']
            total_volume = volume + self.pipe_dictionary['total_volume']
            # total_outer_surface_area = outer_surface_area + self.pipe_dictionary['total_outer_surface_area']

            new_segment = {
                    name : {
                    'material': material,
                    'length': length,
                    'outer_diameter': outer_diameter,
                    'inner_diameter': inner_diameter,
                    'thickness': thickness,
                    'diffusion_path_length': diffusion_path_length,
                    'volume': volume,
                    'outer_surface_area': outer_surface_area,
                    'inner_surface_area': inner_surface_area,

                    },
                    }

            segment_dict = self.pipe_dictionary['segments']
            segment_dict.update(new_segment)
            pipe_dictionary = {
                    'number_segments': self.count,
                    'segment_list': self.segment_list,
                    'total_length':total_length,
                    'total_volume':total_volume,
                    # 'total_outer_surface_area':total_outer_surface_area,
                    'segments': segment_dict
                ,
            }
        else:
            
            pipe_dictionary = {
                    'number_segments': self.count,
                    'segment_list': self.segment_list,
                    'total_length':length,
                    'total_volume':volume,
                    # 'total_outer_surface_area':outer_surface_area,

                'segments': {
                    name : {
                    'material': material,
                    'length': length,
                    'outer_diameter': outer_diameter,
                    'inner_diameter': inner_diameter,
                    'thickness': thickness,
                    'diffusion_path_length': diffusion_path_length,
                    'volume': volume,
                    'outer_surface_area': outer_surface_area,
                    'inner_surface_area': inner_surface_area,
                    },
                },
            }
            
        self.pipe_dictionary = pipe_dictionary
        
        self._calculate_pipe_K_D(pipe_material=material, segment_name = name)

    def _calculate_peak_allowable_gw_concentration_per_segment(self, 
                                    pipe_segment=None,
                                    stagnation_time_hours = None,  
                                    ):
        '''
        Calculates the peak (maximum) concentration in groundwater water for a 
        given a stagnation period, that would not result in a peak concentration 
        in drinking water exceeding the drinking water norm for a single pipe segment.
        Stagnation period default of 8 hours. Peak concentrations in groundwater 
        water and soil added to the pipe_permeability_dict. If the distribution 
        coefficient it unknown for a given chemical, no soil concentration is 
        calculated.
        
        Parameters
        ----------
        stagnation_time_hours: float
            time in hours, default 8 hours
        pipe_segment: string,
            Name of the pipe segment for which the concentrations are calculated
        '''

        drinking_water_norm = self.pipe_permeability_dict['Drinking_water_norm']
        stagnation_time = stagnation_time_hours / 24 # days
        segment_volume = self.pipe_dictionary['segments'][pipe_segment]['volume']
        segment_outer_surface_area = self.pipe_dictionary['segments'][pipe_segment]['outer_surface_area']
        segment_diffusion_path_length = self.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 

        #Risk limit value groundwater
        flux_max_stagnation = ( drinking_water_norm / 1000 * segment_volume /
                        stagnation_time)
        flux_max_stagnation_per_m2 = flux_max_stagnation / segment_outer_surface_area
        
        stagnation_factor = 10 ** max((((self.pipe_permeability_dict['segments'][pipe_segment]['log_Dp'] + 12.5) / 2 + 
                                self.pipe_permeability_dict['segments'][pipe_segment]['log_Kpw']) * 0.73611 + 
                                -1.03574 ), 0)

        concentration_peak_without_stagnation = (flux_max_stagnation_per_m2 * 
                                    segment_diffusion_path_length / 
                                    self.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient'] 
                                    * self.assessment_factor_groundwater)


        concentration_peak_after_stagnation = stagnation_factor * concentration_peak_without_stagnation

        #Risk limit value soil, first check if a distribution coefficient is known
        if math.isnan(self.pipe_permeability_dict['log_distribution_coefficient']):
            concentration_peak_soil = 'log_distribution_coefficient (Kd) unknown'
        else:
            concentration_peak_soil = (10 ** self.pipe_permeability_dict['log_distribution_coefficient'] * 
                                        concentration_peak_after_stagnation * 
                                        self.assessment_factor_soil / self.assessment_factor_groundwater)

        self.pipe_permeability_dict['segments'][pipe_segment]['stagnation_time_hours'] = stagnation_time_hours
        self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_stagnation'] = flux_max_stagnation
        self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_stagnation_per_m2'] = flux_max_stagnation_per_m2
        self.pipe_permeability_dict['segments'][pipe_segment]['stagnation_factor'] = stagnation_factor
        self.pipe_permeability_dict['segments'][pipe_segment]['concentration_peak_without_stagnation'] = concentration_peak_without_stagnation
        self.pipe_permeability_dict['segments'][pipe_segment]['concentration_peak_after_stagnation'] = concentration_peak_after_stagnation
        self.pipe_permeability_dict['segments'][pipe_segment]['concentration_peak_soil'] = concentration_peak_soil

    def _calculate_mean_allowable_gw_concentration_per_segment(self, 
                                    pipe_segment, 
                                    ):
        '''
        Calculates the mean 24 hour concentration in groundwater which would not 
        result in a drinking water concentration exceeding the drinking water
        norm for a single pipe segment. Mean concentrations in groundwater and soil added to the 
        pipe_permeability_dict. If the distribution coefficient it unknown for 
        a given chemical, no soil concentration is calculated.
        
        Parameters
        ----------
        pipe_segment: string,
            Name of the pipe segment for which the concentrations are calculated
        '''
        # Check if the flow rate has been set, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
                             To set flow rate use .set_flow_rate()')
        else: 
            drinking_water_norm = self.pipe_permeability_dict['Drinking_water_norm']
            segment_volume = self.pipe_dictionary['segments'][pipe_segment]['volume']
            segment_outer_surface_area = self.pipe_dictionary['segments'][pipe_segment]['outer_surface_area']
            segment_diffusion_path_length = self.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 

            # 24 hour max flux
            flux_max_per_day = drinking_water_norm / 1000 * self.flow_rate
            flux_max_per_day_per_m2 = flux_max_per_day / segment_outer_surface_area

            concentration_mean = (flux_max_per_day_per_m2 * segment_diffusion_path_length / 
                                        self.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient'] * 
                                        self.assessment_factor_groundwater + drinking_water_norm / 1000)
            
            #Risk limit value soil, first check if there is a distribution coefficient known
            if math.isnan(self.pipe_permeability_dict['log_distribution_coefficient']):
                concentration_mean_soil = 'log_distribution_coefficient (Kd) unknown'
            else:
                concentration_mean_soil = (10 ** self.pipe_permeability_dict['log_distribution_coefficient'] * 
                                            concentration_mean * 
                                            self.assessment_factor_soil / self.assessment_factor_groundwater)

            self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_per_day'] = flux_max_per_day
            self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_per_day_per_m2'] = flux_max_per_day_per_m2
            self.pipe_permeability_dict['segments'][pipe_segment]['concentration_mean'] = concentration_mean
            self.pipe_permeability_dict['segments'][pipe_segment]['concentration_mean_soil'] = concentration_mean_soil

    def _calculate_peak_dw_concentration_per_segment(self, 
                                         pipe_segment=None,
                                        stagnation_time_hours = 8, ):
        '''
        Calculates the peak (maximum) concentration in drinking water for a 
        given a stagnation period given a groundwater concentration, for each pipe segment.
        Stagnation period default of 8 hours. Peak concentrations in drinking 
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
        segment_length = self.pipe_dictionary['segments'][pipe_segment]['length']
        segment_outer_surface_area = self.pipe_dictionary['segments'][pipe_segment]['outer_surface_area']
        segment_diffusion_path_length = self.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length'] 
        concentration_groundwater = self.pipe_permeability_dict['concentration_groundwater']
        segment_diffusion_path_length = self.pipe_dictionary['segments'][pipe_segment]['diffusion_path_length']
        inner_outer_diameter = self.pipe_dictionary['segments'][pipe_segment]['outer_diameter'] - self.pipe_dictionary['segments'][pipe_segment]['thickness'] 
        permeation_coefficient = self.pipe_permeability_dict['segments'][pipe_segment]['permeation_coefficient']
        
        # From equation 4-10 KWR 2016.056
        # LEFT OFF HERE
        concentratie_drinkwater = (permeation_coefficient * 2 * concentration_groundwater * stagnation_time) / (segment_diffusion_path_length * inner_outer_diameter/2)


        flux_max_stagnation = ( concentration_groundwater * segment_volume /
                        stagnation_time)
        
        flux_max_stagnation_per_m2 = flux_max_stagnation / segment_outer_surface_area

        self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_per_day'] = flux_max_stagnation
        self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_per_day_per_m2'] = flux_max_stagnation_per_m2


    def calculate_peak_dw_concentration(self, 
                                        stagnation_time_hours = 8, ):
        '''
        Calculates the peak (maximum) concentration in drinking water for a 
        given a stagnation period given a groundwater concentration.
        Stagnation period default of 8 hours. Peak concentrations in drinking 
        water added to the pipe_permeability_dict.
        
        Parameters
        ----------
        stagnation_time_hours: float
            time in hours, default 8 hours

        '''
        stagnation_time = stagnation_time_hours / 24 # days

        sum_fluxes = 0

        for pipe_segment in self.pipe_dictionary['segment_list']:
            self._calculate_peak_dw_concentration_per_segment(pipe_segment=pipe_segment,
                                    stagnation_time_hours = stagnation_time_hours,  
                                    )
            sum_fluxes += self.pipe_permeability_dict['segments'][pipe_segment]['flux_max_per_day_per_m2']
        
        concentration_pipe_drinking_water = (sum_fluxes / 
                                             self.pipe_dictionary['total_volume'])
        
        self.pipe_permeability_dict['peak_concentration_pipe_drinking_water'] = concentration_pipe_drinking_water

    # AH_todo FUNCTIONS COMPLETE UNTIL HERE
    # LEFT OFF HERE, need to discuss how to loop and calculate the groundwater 
    # concentration on multiple sections

    def calculate_mean_allowable_gw_concentration(self, 
                                    ):
        '''
        Calculates the mean 24 hour concentration in groundwater which would not 
        result in a drinking water concentration exceeding the drinking water
        norm. Mean concentrations in groundwater and soil added to the 
        pipe_permeability_dict. If the distribution coefficient it unknown for 
        a given chemical, no soil concentration is calculated.
        
        '''
        #@MartinvdS, need to adjust this to account for mutliple segments, 
        # right now it doesn't do what its supposed to...

        for pipe_segment in self.pipe_dictionary['segment_list']:
            self._calculate_mean_allowable_gw_concentration_per_segment(pipe_segment=pipe_segment,)


    def calculate_peak_allowable_gw_concentration(self, 
                                    stagnation_time_hours = 8,  
                                    ):
        '''
        Calculates the peak (maximum) concentration in groundwater water for a 
        given a stagnation period that would not result in a peak concentration 
        in drinking water exceeding the drinking water norm for each pipe segment.
        Stagnation period default of 8 hours. Peak concentrations in groundwater 
        water and soil added to the pipe_permeability_dict. If the distribution 
        coefficient it unknown for a given chemical, no soil concentration is 
        calculated.
        
        Parameters
        ----------
        stagnation_time_hours: float
            time in hours, default 8 hours
        pipe_segment: string,
            Name of the pipe segment for which the concentrations are calculated

        '''
        #@MartinvdS, need to adjust this to account for mutliple segments, 
        # right now it doesn't do what its supposed to...

        for pipe_segment in self.pipe_dictionary['segment_list']:
            self._calculate_peak_allowable_gw_concentration_per_segment(pipe_segment=pipe_segment,
                                    stagnation_time_hours = stagnation_time_hours,  
                                    )



    def __str__(self,):
        ''' or override the "print" function '''
        return str(self.pipe_dictionary)
    
    def _set_multidiffusion(self,):
        ''' Function to return diffusion coefficient for multicomponent diffusion'''
