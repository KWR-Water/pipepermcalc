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
    Pipe object class to make segments of the pipe

    Return
    -------
    ....

    '''

    count = 0 # count of pipe segments

    def __init__(self, 
                # some_variable: int or float = None,
                
                ):
        '''
        Parameters
        ----------


        Returns
        -------


        '''       
        segment_list = [] # list of pipe segment names
        self.segment_list = segment_list
    # Reference dictionary with values corresponding to the chemical group numbers
        # Chemical group numbers
        # Expert opinion
        # PAK	1
        # MAK	1
        # ClArom	1
        # ClAlk	1
        # Arom	1
        # Alk	1
        # PCB	2
        # overige	3
        # onbekend	3
        # O2	3
        # Cl	3
        # BDE	3
    
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

    def add_segment(self,
                    name=None,
                    material=None,
                    length=None,
                    diameter=None,
                    thickness=None,):
        ''' Add a segment to the pipe
        
        Parameters
        ----------
        name: string
            name of the pipe segment
        material: enum?? @Bram -> set choice of materials
            e.g. PE40, PE80, PVC, EPDM, rubber etc.
        length: float
            in m ? @AH check units
        diameter: float
            in m ? @AH check units
        thickness: float
            in m ? @AH check units
        
        Returns
        -------
        pipe_dictionary: dictionary
        '''
        xx = 1
        self.count += 1 #count the number of segments created?
        self.segment_list.append(name) 

        if self.count >1:
            volume = math.pi * (diameter / 2) ** 2 * length
            surface_area =  (2 * math.pi * (diameter / 2) * length)

            total_length = length + self.pipe_dictionary['total_length']
            total_volume = volume + self.pipe_dictionary['total_volume']
            total_surface_area = surface_area + self.pipe_dictionary['total_surface_area']

            new_segment = {
                    name : {
                    'material': material,
                    'length': length,
                    'diameter': diameter,
                    'thickness': thickness,
                    },
                    }

            segment_dict = self.pipe_dictionary['segments']
            segment_dict.update(new_segment)
            pipe_dictionary = {
                    'number_segments': self.count,
                    'segment_list': self.segment_list,
                    'total_length':total_length,
                    'total_volume':total_volume,
                    'total_inner_surface_area':total_surface_area,
                    'segments': segment_dict
                ,
            }
        else:
            volume = math.pi * (diameter / 2) ** 2 * length
            surface_area =  (2 * math.pi * (diameter / 2) * length)
            
            pipe_dictionary = {
                    'number_segments': self.count,
                    'segment_list': self.segment_list,
                    'total_length':length,
                    'total_volume':volume,
                    'total_surface_area':surface_area,
                'segments': {
                    name : {
                    'material': material,
                    'length': length,
                    'diameter': diameter,
                    'thickness': thickness,
                    },
                },
            }
        self.pipe_dictionary = pipe_dictionary

# LEFT OFF HERE FINISHED THE DICTIONARIES, 
# NEED TO GO OVER STRUCTURE WITH MARTIN/BRAM/LENNART

    def set_groundwater_conditions(self,
                                   chemical_name=None,                                    
                                   concentration_groundwater=None,
                                   temperature_groundwater=None):
        ''' specify the chemical, concentration and temperature in the groundwater'''
        
        self.concentration_groundwater = concentration_groundwater
        self.temperature_groundwater = temperature_groundwater
        self.chemical_name = chemical_name
        # AH_todo something here to indicate the groundwater conditions have been set?

    def phase_distribute(self,):
        ''' Something to distribute the concentration of the chemical of interest
        to the gas/water/solid phases?'''

    def fetch_chemical_database(self,
                                chemical_name=None,):
        ''' Fetch the pipe and chemical information corresponding to the given pipe 
        material and chemical choice 
        
        Parameters
        ----------
        chemical_name: string 
            @Bram to be replaced by the function restricting input of names

        Returns
        -------
        pipe_permeability_dict: dictionary
            Dictionary of chemical information (solubility, Log Kow, Mwt etc)
        '''
        ppc_database = read_csv(file_path.parent / 'database' / 'ppc_database.csv',  skiprows=[1, 2] ) 

        df = ppc_database[ppc_database['Chemical_name'].str.contains(chemical_name)]

        chemical_dict = df.to_dict('records')[0]

        return chemical_dict

    def calculate_reference_value(self,
                                  reference_pipe_material_dict=None, 
                                  pipe_permeability_dict=None,
                                coefficient_name=None, 
                                coefficient_symbol=None,
                                pipe_material=None,):
        '''Calculate reference log coefficient plastic-water (log pw) 

        Parameters
        ----------
        reference_pipe_material_dict: dictionary
            Dictionary with reference values for different pipe materials
            corresponding to the chemical group numbers
        pipe_permeability_dict: dictionary
        coefficient_name: string 
        coefficient_symbol: string
        pipe_material: string
            Choice of pipe material: PE40, PE80, SBR, EPDM

        Returns
        -------
        log_pw_ref: float
            Reference log coefficient for given plastic/plastic-water
        '''

        a_ref = reference_pipe_material_dict[pipe_material]['ref_log_'+coefficient_symbol+ '_a'][pipe_permeability_dict['chemical_group_number']]
        b_ref = reference_pipe_material_dict[pipe_material]['ref_log_'+coefficient_symbol+ '_b'][pipe_permeability_dict['chemical_group_number']]
        log_pw_ref = a_ref * pipe_permeability_dict[coefficient_name] + b_ref
        return log_pw_ref

    def correct_for_temperature(self,
                                pipe_permeability_dict=None, 
                        temperature_groundwater=None,
                        coefficient_name=None, 
                        a_dh=None,
                        b_dh=None,
                        ):
        '''Temperature correction for Log pw, 
        From ppc database, sheet Material

        Parameters
        ----------


        Returns
        -------

        '''
        R = 0.008314 #universal gas constant
        reference_temperature = 25 # deg. C
        dh = a_dh * math.log10(pipe_permeability_dict[coefficient_name]) + b_dh
        f_temp = dh / (R * math.log(10)) * (1 / (reference_temperature + 273) - 1 / (temperature_groundwater + 273))
        return f_temp

    def other_correction(self,
                         pipe_permeability_dict=None, 
                        a_c=None,
                        Cref_Sw=None):
        '''Correction factor, @ah rename this function when checked what it is

        Parameters
        ----------


        Returns
        -------
        '''

        # from ppc_database material K27-29
        # a_c = 0.103965019849463
        # Cref_Sw = 1.000
        Cg_Sw = pipe_permeability_dict['concentration_groundwater'] / pipe_permeability_dict['solubility']
        f_conc = a_c * (Cg_Sw - Cref_Sw)
        return f_conc

    def correct_for_age(self,):
        '''Age correction @ah check this with @MartinvdS
        for now a placeholder function if we ever wanted to add this later '''
        f_age = 0.000
        return f_age

    def calculate_pipe_K_D(self,
                        chemical_name=None, 
                        pipe_material=None, ):
        ''' Fetch the pipe and chemical information corresponding to the given pipe 
        material and chemical choice 

        Parameters
        ----------
        chemical_name: string 
            @Bram to be replaced by the function restricting input of names
        pipe_material: string
            Choice of pipe material: PE40, PE80, SBR, EPDM
        temperature_groundwater: float
            Temperature of groundwater, in degrees Celcius
        concentration_groundwater: float
            concentration of given chemical in groundwater, in mg/L

        Returns
        -------
        '''
        
        reference_pipe_material_dict = self.reference_pipe_material_dict

        pipe_permeability_dict = self.fetch_chemical_database(chemical_name=chemical_name)

        #AH_todo add a check here to see if groundwater conditions have been set
        pipe_permeability_dict['concentration_groundwater'] = self.concentration_groundwater
        pipe_permeability_dict['temperature_groundwater'] = self.temperature_groundwater

        # calculate reference log K plastic-water (log kpw) 
        log_Kpw_ref = self.calculate_reference_value(reference_pipe_material_dict=self.reference_pipe_material_dict, 
                            pipe_permeability_dict=pipe_permeability_dict,
                            coefficient_name ='octanol_water_partitioning_coefficient', 
                            coefficient_symbol ='K',
                            pipe_material=pipe_material)
        pipe_permeability_dict['log_Kpw_ref'] = log_Kpw_ref

        f_Ktemp = self.correct_for_temperature(pipe_permeability_dict=pipe_permeability_dict, 
                        temperature_groundwater=self.temperature_groundwater, 
                        coefficient_name = 'solubility',
                            a_dh = 7.92169801506708,
                            b_dh = -17.1875608983359,
                        )

        f_Kconc = self.other_correction(pipe_permeability_dict=pipe_permeability_dict,
                                a_c = 0.103965019849463,
                                Cref_Sw = 1.000)
        
        f_Kage = self.correct_for_age()

        log_Kpw = log_Kpw_ref + f_Ktemp + f_Kconc + f_Kage

        pipe_permeability_dict['f_Ktemp'] = f_Ktemp   
        pipe_permeability_dict['f_Kconc'] = f_Kconc
        pipe_permeability_dict['log_Kpw'] = log_Kpw


        # calculate reference log D plastic-water (log Dpw) 
        log_Dp_ref= self.calculate_reference_value(reference_pipe_material_dict=self.reference_pipe_material_dict, 
                            pipe_permeability_dict=pipe_permeability_dict,
                            coefficient_name ='molecular_weight', 
                            coefficient_symbol ='D',
                            pipe_material=pipe_material)
        pipe_permeability_dict['log_Dp_ref'] = log_Dp_ref
        
        f_Dtemp = self.correct_for_temperature(pipe_permeability_dict=pipe_permeability_dict, 
                        temperature_groundwater=self.temperature_groundwater, 
                        coefficient_name ='molecular_weight', 
                            a_dh = 61.8565740136974,
                            b_dh = -78.9191401984509,
                        )

        f_Dconc = self.other_correction(pipe_permeability_dict=pipe_permeability_dict,
                                a_c = 0.784077209735583,
                                Cref_Sw = 0.5)
        
        f_Dage = self.correct_for_age()

        log_Dp = log_Dp_ref + f_Dtemp + f_Dconc + f_Dage

        pipe_permeability_dict['f_Dtemp'] = f_Dtemp    
        pipe_permeability_dict['f_Dconc'] = f_Dconc
        pipe_permeability_dict['log_Dp'] = log_Dp

        pipe_permeability_dict['Ppw'] = 24 * 60 * 60 * (10 ** log_Dp) * 10 ** log_Kpw

        self.pipe_permeability_dict = pipe_permeability_dict
        return pipe_permeability_dict #@AH necessary to return this?

    def calculate_max_dw_concentration(self,):
        ''' Calculate the peak/max concentration in drinking water'''
        
    def print_pipe_segment_information(self,):
        ''' or override the "print" function '''

    def print_chemical_pipe_information(self,):
        ''' Print informaiton about the chemical/pipe
        e.g. K, D etc '''

    def _set_diffusion_coefficient(self,):
        ''' Functiont to override the diffusion coefficient from the database, 
        private function...'''
    
    def _set_multidiffusion(self,):
        ''' Function to return diffusion coefficient for multicomponent diffusion'''

