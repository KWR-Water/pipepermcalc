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

from project_path import file_path, module_path

from pipepermcalc.segment import * 

class Pipe:
    '''
    Pipe object class to make segments of the pipe and calculate the peak and
    mean concentration of a chemical in groundwater and soil.

    Attributes
    -------
    segment_list: list
        list of the pipe segment objects which make up the pipe
   
    pipe_dictionary: dictionary
        number_segments: float
            total number of pipe segments
        segment_list: list
            list of the names of the different pipe sements
        total_length: float
            sum of the lengths of the pipe segments
        total_volume: float
            sum of the volumes of the pipe segments
        total_outer_surface_area: float,
            sum of the surface area of the pipe segments
        flow_rate: float
            flow rate in the pipe, m3/day 
        segments: dictionary
            dictionary of the individual pipe segments, containing the 
            segment material, segment length (m), inner_diameter (m), thickness (m), 
            volume (m3) and surface area (m2), permeation_direction and
            diffusion_path_length (m).

    pipe_permeability_dict: dictionary
        CAS_number: string
            CAS is a unique identification number assigned by the Chemical 
            Abstracts Service (CAS)
        chemical_name_NL: string
            Name of the chemical given in Dutch
        chemical_name: string
            @Bram to be replaced by the function restricting input of names
        molecular_weight: float
            Mass of one mole of a given chemical, g/mol
        solubility:	float
            solubility of given chemical in water, g/m3
        log_octanol_water_partitioning_coefficient:	float,
            Partition coefficient for the two-phase system consisting of 
            n-octanol and water, Log Kow, [-]
        log_distribution_coefficient: float
            Ratio of the amount of chemical  adsorbed onto soil per amount 
            of water, m3/g
        chemical_group: string
            Grouping of chemicals (expert opinion) with similar properties
            for permeation: Group 1: PAK, MAK, ClArom, ClAlk, Arom, Alk, 
            Group 2: PCB, Group 3: overig, onbekend, O2, Cl, BDE. See KWR 2016.056
        chemical_group_number: integer
            Integer corresponding to the chemical group 
        molecular_volume: float
            Volume occupied by one mole of the substance at a given 
            temperature and pressure, cm3/mol
        Drinking_water_norm: float
            Concentration allowable in the Dutch Drinking water decree, g/m3
        concentration_groundwater: float
            Concentration of the given chemical in groundwater, g/m3
        temperature_groundwater: float
            Temperature of the groundwater, degrees Celcius
        
        tolerance: float 
            the allowable difference between the calculated and actual drinking water concentration
        relaxatoin_factor: float
            used to iterate and calculate the new drinking water concentration
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme
                    
        #ah_todo update these after finished these functions.
            
        # stagnation_time_hours:float
        #     Time in hours which water in pipe is stagnant, hours.
        # flux_max_stagnation: float
        #     Maximum flux during stagnation period to remain below the drinking 
        #     water standard, mg/day
        # flux_max_stagnation_per_m2: float
        #     Maximum flux per square meter (surface area) of pipe during 
        #     stagnation period to remain below the drinking water standard, g/day
        # flux_max_per_day: float
        #     Maximum flux in a day (24 hours) to remain below the drinking 
        #     water standard, g/day
        # flux_max_per_day_per_m2: float
        #     Maximum flux per square meter (surface area) of pipe in one day 
        #     (24 hours) to remain below the drinking water standard, g/day        
        
        # concentration_gw_peak_without_stagnation: float
        #     Concentration in groundwater which, wihtout a stagnation period, 
        #     would not result in a peak concentration in drinking water exceeding 
        #     the drinking water norm, g/m3
        # concentration_gw_peak_after_stagnation: float
        #     Concentration in groundwater which, after a stagnation period, 
        #     would not result in a peak concentration in drinking water exceeding 
        #     the drinking water norm, g/m3
        # concentration_peak_soil: float
        #     Concentration in soil which, after a stagnation period, 
        #     would not result in a peak concentration in drinking water exceeding 
        #     the drinking water norm, mg/kg
        # concentration_gw_mean:float
        #     Mean concentration in groundwater which would would not result in 
        #     a mean daily (24 horus) concentration in drinking water exceeding 
        #     the drinking water norm, g/m3
        # concentration_mean_soil: float     
        #     Mean concentration in soil which would would not result in 
        #     a mean daily (24 horus) concentration in drinking water exceeding 
        #     the drinking water norm, mg/kg
    '''
    #ah_todo change input variables to restrict the type (e.g. only float, 
    # only integer, only positive values etc)
    
    def __init__(self, 
                 segment_list,
                ):
        '''
        segment_list: list
            list of the segments objects
            
        '''       
        self.segment_list = segment_list
        self._groundwater_conditions_set = False
        self._flow_rate_set = False

        sum_total_volume = 0
        sum_total_length = 0
        for segment in segment_list:
            sum_total_length += segment.length
            sum_total_volume += segment.volume

        self.total_length = sum_total_length
        self.total_volume = sum_total_volume


    def check_input_values(self, check_values):
        '''Checks that the input values are > 0
        
        Parameters
        ----------
        check_values: list,
            List (string) of the variables to check 
        '''
        
        for check_value in check_values:
            value = getattr(self, check_value)
            if value <= 0:
                raise ValueError(f'{check_value} must be > 0 ')
        

    def set_groundwater_conditions(self,
                                   chemical_name=None,                                    
                                   concentration_groundwater=None,
                                   temperature_groundwater=None):
        ''' 
        Specifies the chemical of interest, concentration and temperature in the 
        groundwater and returns the parameters as attributes of the class. 
        Calculates the segment permeation parameters based on the groundwater 
        conditions.
        
        Parameters
        ----------
        chemical_name: string
            @Bram to be replaced by the function restricting input of names
        concentration_groundwater: float
            Concentration of the given chemical in groundwater, g/m3
        temperature_groundwater: float
            Temperature of the groundwater, degrees Celcius
        '''


        self.concentration_groundwater = concentration_groundwater
        self.temperature_groundwater = temperature_groundwater
        # Checks here that input concentration and temperature > 0
        check_values = ['concentration_groundwater', 'temperature_groundwater',]
        self.check_input_values(check_values)

        self.chemical_name = chemical_name
        self._groundwater_conditions_set = True

        self.pipe_permeability_dict = self._fetch_chemical_database(chemical_name=self.chemical_name)
        self.pipe_permeability_dict['chemical_name'] = self.chemical_name
        self.pipe_permeability_dict['concentration_groundwater'] = self.concentration_groundwater
        self.pipe_permeability_dict['temperature_groundwater'] = self.temperature_groundwater

        for segment in self.segment_list:
            segment._calculate_pipe_K_D(pipe_permeability_dict=self.pipe_permeability_dict, 
                                 _groundwater_conditions_set=self._groundwater_conditions_set, )

    

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

        # Checks here that input flow_rate > 0
        check_values = ['flow_rate', ]
        self.check_input_values(check_values)

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
        
        ppc_database = read_csv(module_path / 'database' / 'ppc_database.csv',  skiprows=[1, 2] ) 

        df = ppc_database[ppc_database['chemical_name'].str.contains(chemical_name)]

        chemical_dict = df.to_dict('records')[0]
        # convert to g/m3
        chemical_dict['Drinking_water_norm'] = chemical_dict['Drinking_water_norm']/1000 

        return chemical_dict


    def calculate_mean_dw_concentration(self, 
                                        tolerance = 0.01, #ah_todo should we have these as defaults?
                                        relaxation_factor = 0.5,
                                        max_iterations = 1000):
        '''
        Calculates the mean concentration in drinking water for a 24 hour period
        given a groundwater concentration. Mean concentrations in drinking 
        water added to the pipe_permeability_dict.
        
        Parameters
        ----------
        tolerance: float 
            the allowable difference between the calculated and actual drinking water concentration
        relaxatoin_factor: float
            used to iterate and calculate the new drinking water concentration
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        '''
        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor
        
        # Checks here that input 'max_iterations', 'tolerance', 'relaxation_factor' > 0
        check_values = ['max_iterations', 'tolerance', 'relaxation_factor' ]
        self.check_input_values(check_values)


        # Check if the flow rate has been set, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')
        else: 
            concentration_drinking_water = 0 #initial guess for drinking water 
            counter = 0

            while True:    

                sum_mass_segment = 0

                for segment in self.segment_list:
                    segment._calculate_mean_dw_mass_per_segment(pipe_permeability_dict=self.pipe_permeability_dict,
                                                                concentration_drinking_water = concentration_drinking_water,
                                                                _groundwater_conditions_set = self._groundwater_conditions_set,
                                                                flow_rate = self.flow_rate)

                    sum_mass_segment += segment.mass_chemical_drinkwater
            
                concentration_pipe_drinking_water = (sum_mass_segment / 
                                                self.flow_rate) #volume of water consumed in 1 day = flow rate
            
                
                counter +=1
                
                if abs(1 - concentration_drinking_water / concentration_pipe_drinking_water) / relaxation_factor <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                else:
                    concentration_drinking_water = (relaxation_factor 
                                                    * concentration_pipe_drinking_water 
                                                    + (1 - relaxation_factor) 
                                                    * concentration_drinking_water)
                # if counter % 10 ==0 : print(concentration_drinking_water) #for debugging
                
            self.pipe_permeability_dict['mean_concentration_pipe_drinking_water'] = concentration_pipe_drinking_water


    def calculate_peak_dw_concentration(self, 
                                        stagnation_time_hours = 8, 
                                        tolerance = 0.01, #ah_todo should we not have these as defaults?
                                        relaxation_factor = 0.5,
                                        max_iterations = 1000):

        '''
        Calculates the peak (maximum) concentration in drinking water for a 
        given a stagnation period given a groundwater concentration.
        Stagnation period default of 8 hours. Peak concentrations in drinking 
        water added to the pipe_permeability_dict.
        
        Parameters
        ----------
        stagnation_time_hours: float
            time in hours, default 8 hours
        tolerance: float 
            the allowable difference between the calculated and actual drinking water concentration
        relaxatoin_factor: float
            used to iterate and calculate the new drinking water concentration
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        '''

        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor
        self.stagnation_time_hours = stagnation_time_hours

        # Checks here that input 'stagnation_time_hours', 'max_iterations', 
        # 'tolerance', 'relaxation_factor' > 0
        check_values = ['stagnation_time_hours', 'max_iterations', 'tolerance', 'relaxation_factor' ]
        self.check_input_values(check_values)

        # Check if the flow rate has been set, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')
        else: 
            concentration_drinking_water = 0.0 #initial guess for drinking water
            counter = 0

            while True:    

                sum_mass_segment = 0

                for segment in self.segment_list:
                    segment._calculate_peak_dw_mass_per_segment(pipe_permeability_dict=self.pipe_permeability_dict,
                                                                concentration_drinking_water = concentration_drinking_water,
                                                                _groundwater_conditions_set = self._groundwater_conditions_set,
                                                                stagnation_time_hours = stagnation_time_hours, 
                                                                flow_rate = self.flow_rate)

                    sum_mass_segment += segment.mass_chemical_drinkwater
            
                concentration_pipe_drinking_water = (sum_mass_segment / 
                                                self.total_volume ) #volume of water in the pipe during stagnation time
                
                counter +=1
                
                if abs(1 - concentration_drinking_water / concentration_pipe_drinking_water) / relaxation_factor <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                else:
                    concentration_drinking_water = (relaxation_factor 
                                                    * concentration_pipe_drinking_water 
                                                    + (1 - relaxation_factor) 
                                                    * concentration_drinking_water)
                # if counter % 100 ==0 : print(concentration_drinking_water) #for debugging
                
            self.pipe_permeability_dict['peak_concentration_pipe_drinking_water'] = concentration_pipe_drinking_water


    def calculate_mean_allowable_gw_concentration(self, #ah_todo write test
                                        concentration_drinking_water,
                                        chemical_name,
                                        temperature_groundwater,
                                        tolerance = 0.01, #ah_todo should we have these as defaults?
                                        relaxation_factor = 0.5,
                                        max_iterations = 1000
                                        ):
        '''
        Calculates the mean 24 hour concentration in groundwater which would not 
        result in a drinking water concentration exceeding the drinking water
        norm. Mean concentrations in groundwater and soil added to the 
        pipe_permeability_dict. If the distribution coefficient it unknown for 
        a given chemical, no soil concentration is calculated.

        Parameters
        ----------
        concentration_drinking_water: float
            Concentration in the drinking water for which to calculate the mean 
            allowable groundwater concentration, g/m^3
        tolerance: float 
            the allowable difference between the calculated and actual drinking water concentration
        relaxatoin_factor: float
            used to iterate and calculate the new drinking water concentration
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme
        
        '''
        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor
        
        # Checks here that input 'max_iterations', 'tolerance', 'relaxation_factor' > 0
        check_values = ['max_iterations', 'tolerance', 'relaxation_factor' ]
        self.check_input_values(check_values)


        # Check if the flow rate has been set, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')
        else: 
            pipe_permeability_dict = self._fetch_chemical_database(
                                            chemical_name=chemical_name)

            # calculate initial guess for gw concentration
            sum_KDA_d = 0
            for segment in self.segment_list:
                # calculate the sum of the Kpw * DP * SA / d for all pipe segments
                log_Dp_ref = segment._calculate_ref_logD(pipe_permeability_dict=pipe_permeability_dict, )
                log_Kpw_ref = segment._calculate_ref_logK(pipe_permeability_dict=pipe_permeability_dict, )
                
                Dp_ref = 10 ** log_Dp_ref
                Kpw_ref = 10 ** log_Kpw_ref

                permeation_coefficient_ref = (24 * 60 * 60 * Dp_ref * Kpw_ref)

                sum_KDA_d_segment = (permeation_coefficient_ref * segment.permeation_surface_area 
                                    / segment.diffusion_path_length )

                sum_KDA_d += sum_KDA_d_segment

            concentration_groundwater = (concentration_drinking_water * (1
                                         + self.flow_rate * segment.assessment_factor_groundwater ) 
                                            / sum_KDA_d )
            counter = 0

            while True:
                self.set_groundwater_conditions(chemical_name=chemical_name, 
                                            temperature_groundwater=temperature_groundwater, 
                                            concentration_groundwater=concentration_groundwater, 
                                            )
                sum_mass_segment = 0

                # mass of chemical in pipe water to meet drinking water norm
                mass_drinkingwater_norm = concentration_drinking_water * self.flow_rate

                for segment in self.segment_list:
                    segment._calculate_mean_dw_mass_per_segment(pipe_permeability_dict=self.pipe_permeability_dict,
                                                                    concentration_drinking_water = concentration_drinking_water,
                                                                    _groundwater_conditions_set = self._groundwater_conditions_set,
                                                                    flow_rate = self.flow_rate)
                    sum_mass_segment += segment.mass_chemical_drinkwater

                counter +=1

                if abs(1 - mass_drinkingwater_norm / sum_mass_segment) <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                else:
                    new_groundwater = concentration_groundwater * (1 - relaxation_factor + relaxation_factor * (mass_drinkingwater_norm / sum_mass_segment))
                    concentration_groundwater = new_groundwater
                    # if counter % 100 ==0 : print(concentration_drinking_water) #for debugging



    def calculate_peak_allowable_gw_concentration(self, #ah_todo write test
                                    concentration_drinking_water,
                                    chemical_name,
                                    temperature_groundwater,
                                    stagnation_time_hours = 8,
                                    tolerance = 0.01, #ah_todo should we have these as defaults?
                                    relaxation_factor = 0.5,
                                    max_iterations = 1000

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
        tolerance: float 
            the allowable difference between the calculated and actual drinking water concentration
        relaxatoin_factor: float
            used to iterate and calculate the new drinking water concentration
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        '''
        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor
        self.stagnation_time_hours = stagnation_time_hours

        # Checks here that input 'stagnation_time_hours', 'max_iterations', 
        # 'tolerance', 'relaxation_factor' > 0
        check_values = ['stagnation_time_hours', 'max_iterations', 'tolerance', 'relaxation_factor' ]
        self.check_input_values(check_values)

        # Check if the flow rate has been set, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')
        else: 

            pipe_permeability_dict = self._fetch_chemical_database(chemical_name=chemical_name)
            self.stagnation_time = stagnation_time_hours / 24

            # calculate initial guess for gw concentration
            sum_KDA_d = 0
            for segment in self.segment_list:
                # calculate the sum of the Kpw * DP * SA *f_stag / d for all pipe segments
                log_Dp_ref = segment._calculate_ref_logD(pipe_permeability_dict=pipe_permeability_dict, )
                log_Kpw_ref = segment._calculate_ref_logK(pipe_permeability_dict=pipe_permeability_dict, )
                
                Dp_ref = 10 ** log_Dp_ref
                Kpw_ref = 10 ** log_Kpw_ref

                #stagnation factor with reference values for LogDp and LogKpw
                stagnation_factor = 10 ** max((((log_Dp_ref + 12.5) / 2 + 
                                    log_Kpw_ref) * 0.73611 + 
                                    -1.03574 ), 0)            

                permeation_coefficient_ref = (24 * 60 * 60 * Dp_ref * Kpw_ref)

                sum_KDA_d_segment = (permeation_coefficient_ref * segment.permeation_surface_area 
                                    * stagnation_factor
                                    / segment.diffusion_path_length )

                sum_KDA_d += sum_KDA_d_segment

            # initial guess concentration in groundwater
            concentration_groundwater = concentration_drinking_water * (1 
                                        + self.total_volume * segment.assessment_factor_groundwater 
                                        / self.stagnation_time / sum_KDA_d)
            
            counter = 0

            while True:
                self.set_groundwater_conditions(chemical_name=chemical_name, 
                                            temperature_groundwater=temperature_groundwater, 
                                            concentration_groundwater=concentration_groundwater, 
                                            )
                sum_mass_segment = 0

                # mass of chemical in pipe water to meet drinking water norm
                mass_drinkingwater_norm = (concentration_drinking_water * self.total_volume)
                
                for segment in self.segment_list:
                    segment._calculate_peak_dw_mass_per_segment(pipe_permeability_dict=self.pipe_permeability_dict,
                                                                    concentration_drinking_water = concentration_drinking_water,
                                                                    _groundwater_conditions_set = self._groundwater_conditions_set,
                                                                    stagnation_time_hours = stagnation_time_hours,
                                                                    flow_rate = self.flow_rate)
                    sum_mass_segment += segment.mass_chemical_drinkwater

                counter +=1
                
                if abs(1 - mass_drinkingwater_norm / sum_mass_segment) <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                else:
                    new_groundwater = concentration_groundwater * (1 - relaxation_factor + relaxation_factor * (mass_drinkingwater_norm / sum_mass_segment))
                    concentration_groundwater = new_groundwater
                    # if counter % 100 ==0 : print(concentration_drinking_water) #for debugging
                

