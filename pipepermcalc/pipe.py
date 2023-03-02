#%% ----------------------------------------------------------------------------
# A. Hockin, January 2023
# KWR 403230-003
# Pipe permeation calculator
# With Martin vd Schans, Bram Hillebrand, Lennart Brokx
#
# ------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import fuzzywuzzy.fuzz as fwf
import fuzzywuzzy.process as fwp

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
    _conditions_set: Boolean
        Default False, True when the groundwater conditions have been set.
    _flow_rate_set: Boolean
        Default False, True when the flow rate has been set.
    total_volume: float
        Total volume of the pipe, summed from the pipe segments, m3
    total_length': float
        Total length of the pipe, summed from the pipe segments, m
    flow_rate: float
        flow_rate through pipe. Default of 0.5 m3/day.    
    CAS_number: string
        CAS is a unique identification number assigned by the Chemical 
        Abstracts Service (CAS)
    chemical_name_EN: string
        Name of the chemical in English
    chemical_name_NL: string
        Name of the chemical in Dutch
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
        temperature and pressure, cm3/mol.
    Drinking_water_norm: float
        Concentration allowable in the Dutch Drinking water decree, g/m3.
    concentration_groundwater: float
        Concentration of the given chemical in groundwater, g/m3.
    temperature_groundwater: float
        Temperature of the groundwater, degrees Celcius.
    tolerance: float 
        The allowable difference between the calculated and actual drinking 
        water concentration, [-].
    relaxation_factor: float
        Used to iterate and calculate the new drinking water concentration, 
        recommended 0.3-0.7 [-].
    max_iterations: int
        Maximum number of iterations allowed in the optimization scheme.                    
    stagnation_time: float
        Time in seconds which water in pipe is stagnant, unit of seconds. The 
        stagnation factor is only valid for a stagnation time of 8 hours 
        (28800 seconds), therefore using another other stagnation time is not advised.
    concentration_peak_allowable_groundwater: float
        Concentration in groundwater which, after a stagnation period, 
        would not result in a peak concentration in drinking water exceeding 
        the drinking water norm, g/m3.
    concentration_mean_allowable_groundwater: float
        Mean concentration in groundwater which would would not result in 
        a mean daily (24 hours) concentration in drinking water exceeding 
        the drinking water norm, g/m3.
    mean_concentration_pipe_drinking_water: float
        Calculates the mean concentration in drinking water for a 24 hour period
        given a groundwater concentration.
    peak_concentration_pipe_drinking_water: float
        Calculates the peak (maximum) concentration in drinking water for a 
        given a stagnation period given a groundwater concentration.
    
    Note
    ----
    All parameters are in SI units: m, m2, g/m3 (equivalent to mg/L), seconds.


        '''
    #ah_todo change input variables to restrict the type (e.g. only float, 
    # only integer, only positive values etc)
    
    #Constants for iterative calculations
    TOLERANCE_DEFAULT = 0.01
    RELAXATION_FACTOR_DEFAULT = 0.5
    MAX_ITERATIONS_DEFAULT = 1000

    ppc_database = pd.read_csv(module_path / 'database' / 'ppc_database.csv',  skiprows=[1, 2] ) 

    #dictionary to check input parameter value and dtype
    parameter_validation_dictionary = \
        {
        'name':     {'value_dtype': [str]}, 
        'material': {'str_options': ['PE40', 'PE80', 'PVC', 'EPDM'],
                    'value_dtype': [str]}, 
        'permeation_direction': {'str_options': ['perpendicular', 'parallel'],
                    'value_dtype': [str]}, 
        'length': {'min_value': 0, 
                'value_dtype': [float, int]}, 
        'inner_diameter': {'min_value': 0, 
                        'value_dtype': [float, int]}, 
        'wall_thickness': {'min_value': 0, 
                        'value_dtype': [float, int]},  
        'diffusion_path_length': {'min_value': 0, 
                                'value_dtype': [float, int]},
        'max_iterations': {'min_value': 0, 
                        'value_dtype': [float, int]}, 
        'tolerance': {'min_value': 0,
                    'max_value': 1, 
                    'value_dtype': [float]},  
        'relaxation_factor': {'min_value': 0, 
                            'max_value': 1, 
                            'value_dtype': [float]},  
        'stagnation_time': {'min_value': 0, 
                            'value_dtype': [float, int]},  
        'flow_rate': {'min_value': 0, 
                    'value_dtype': [float, int]},  
        'concentration_groundwater': {'min_value': 0, 
                                    'value_dtype': [float, int]},  
        'temperature_groundwater': {'min_value': 0, 
                                    'value_dtype': [float, int]},  
        'concentration_drinking_water': {'min_value': 0, 
                                        'value_dtype': [float, int]},  
        'chemical_name': {'value_dtype': [str]},  
        'language': {'str_options': ['NL', 'EN'],
                    'value_dtype': [str]}, 
        }
    
    def __init__(self, 
                 segment_list,
                ):
        '''
        segment_list: list
            list of the segments objects
            
        ''' 

        self.segment_list = segment_list
        self._conditions_set = False
        self._flow_rate_set = False
        self._is_validated = False

        sum_total_volume = 0
        sum_total_length = 0
        for segment in segment_list:
            sum_total_length += segment.length
            sum_total_volume += segment.volume

        self.total_length = sum_total_length
        self.total_volume = sum_total_volume
    

    def _validate_object(self, check_object):
        ''' Check that the input parameters are valid values and types for the 
        input object'''

        for k, v in self.parameter_validation_dictionary.items():
            # print(getattr(seg1, k))
            if hasattr(check_object, k):
                if type(getattr(check_object, k)) not in v['value_dtype']:
                    raise ValueError(f"Invalid value ~{getattr(check_object, k)}~ for parameter {k}. Input value should be a {v['value_dtype']}.")
                if 'min_value' in v.keys():
                    if getattr(check_object, k) < v['min_value']:
                        raise ValueError(f"Invalid value {getattr(check_object, k)} for parameter {k}. Input value should be a > {v['min_value']}.")
                if 'max_value' in v.keys():
                    if getattr(check_object, k) > v['max_value']:
                        raise ValueError(f"Invalid value ~{getattr(check_object, k)}~ for parameter {k}. Input value should be a < {v['max_value']}.")
                if 'str_options' in v.keys():
                    if getattr(check_object, k) not in v['str_options']:
                        raise ValueError(f"Invalid value ~{getattr(check_object, k)}~ for parameter {k}. Input value should be one of {v['str_options']}.")


    def validate_input_parameters(self,):
        ''' Check that the input parameters are valid values and types for the 
        Pipe and Segment objects'''

        for segment in self.segment_list:
            self._validate_object(segment)
        
        self._validate_object(self)

        self._is_validated=True


    def _fuzzy_min_score(self, chemical_name): #From Vincent Post
        """
        This function calculates the minimum score required for a valid
        match in fuzzywuzzy's extractOne function. The minimum score depends
        on the length of 's' and is calculated based on the string lengths and
        scores in the DEFAULT_MINSCORES dictionary.

        Parameters
        ----------
        chemical_name : str
            String for which the minimum score must be determined, in our case 
            the chemical name

        Returns
        -------
        result : float
            The minimum score for 's'.
        """
        DEFAULT_FUZZY_MINSCORES = {1: 100, 3: 100, 4: 90, 5: 85, 6: 80, 8: 75}

        xp = list(DEFAULT_FUZZY_MINSCORES.keys()) #@martin comment: #@@ maybe it's easier to enter a list directly, instead of extract the dictionary.
        fp = [v for v in DEFAULT_FUZZY_MINSCORES.values()]
        # Use the interp function from NumPy. By default this function
        # yields fp[0] for x < xp[0] and fp[-1] for x > xp[-1]
        return np.interp(len(chemical_name), xp, fp)


    def _extract_matching_chemical_name(self, chemical_name, database):
        ''' Search and extract the highest matching chemical name from the database for the given input
        Parameters
        ----------
        chemical_name : str
            String for which the minimum score/highest matching chemical name 
            must be determined.

        Returns
        -------
        matching_chemical_name : str
            Name with the highest match for the input chemical name from the 
            database.
        '''

        # Exctract the highest scoring chemical name matching the 
        minscore = self._fuzzy_min_score(chemical_name=chemical_name)

        # Return only the highest scoring item
        fuzzy_score = fwp.extractOne(
            query=chemical_name,
            choices=database,
            scorer=fwf.token_sort_ratio,
            score_cutoff=minscore,
        )
        
        matching_chemical_name = fuzzy_score[0]

        return matching_chemical_name


    def _fetch_chemical_database(self,
                                chemical_name=None,
                                suppress_print=False,
                                language = 'NL'
                                ):
        ''' 
        Fetch the pipe and chemical information corresponding to the given 
        pipe material and chemical choice and creates a dictionary 
        pipe_permeability_dict which consists of chemical and permeability 
        related coefficients.

        Parameters
        ----------
        chemical_name: str
            Name of the chemical for which to calculate the permeation
        suppress_print: Boolean
            Suppress printing the chemical name and matching name, e.g. in loop calculations
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible
        '''
        
        database = list(self.ppc_database['chemical_name_'+language])
        
        matching_chemical_name = self._extract_matching_chemical_name(chemical_name=chemical_name, 
                                             database=database)
        
        #ah_todo, @Bram what kind of check here for the chemical name match
        if suppress_print:
            pass
        else:
            print("Input chemical name: ", str(chemical_name), "- Matched chemical name: ", str(matching_chemical_name))

        df = self.ppc_database.loc[self.ppc_database['chemical_name_'+language] == matching_chemical_name]
        pipe_permeability_dict = df.to_dict('records')[0]

        #assign dict items as attribute of class
        for k, v in pipe_permeability_dict.items():
            setattr(self, k, v)


    def set_conditions(self,
                    chemical_name=None,                                    
                    concentration_groundwater=None,
                    temperature_groundwater=None, 
                    flow_rate=None,
                    concentration_drinking_water=None, #ah_todo finish  add to param list
                    stagnation_time = 8 * 60 * 60,
                    suppress_print = False, 
                    language = 'NL'
                    ):
        ''' 
        Specifies the chemical of interest, concentration and temperature in the 
        groundwater and returns the parameters as attributes of the class. 
        Calculates the segment permeation parameters based on the groundwater 
        conditions.
        
        Parameters
        ----------
        chemical_name: string
            Name of the chemical for which to calculate the permeation
        concentration_groundwater: float
            Concentration of the given chemical in groundwater, g/m3
        temperature_groundwater: float
            Temperature of the groundwater, degrees Celcius
        stagnation_time: float
            Time in seconds which water in pipe is stagnant, unit of seconds. The 
            stagnation factor is only valid for a stagnation time of 8 hours 
            (28800 seconds), therefore using another other stagnation time is not advised.
        suppress_print: Boolean
            Suppress printing the chemical name and matching name, e.g. in loop calculations
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible

        '''
        self.chemical_name = chemical_name
        self.concentration_groundwater = concentration_groundwater
        self.temperature_groundwater = temperature_groundwater
        self.stagnation_time = stagnation_time

        #return to these checks...
        self._conditions_set = True

        self.flow_rate = flow_rate
        self._flow_rate_set = True

        self._fetch_chemical_database(chemical_name=self.chemical_name, 
                                        suppress_print=suppress_print, 
                                        language=language)

        # The default value for the concentration_drinking_water is the drinking water norm
        if concentration_drinking_water is None:
            self.concentration_drinking_water = self.Drinking_water_norm
        else: 
            self.concentration_drinking_water = concentration_drinking_water

        # #ah_todo move to the individual calculations 
        # @martin, if we leave this here, then it is ok, we validate after and 
        # any errors are caught before the calculations and then we can override
        # see tag ah_removehere
        for segment in self.segment_list:          
            segment._calculate_pipe_K_D(pipe = self, 
                                        _conditions_set=self._conditions_set, )


    def view_database_chemical_names(self, language='NL'):
        '''function to view a list of the possible chemical names

        Parameters
        ----------
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible
        '''

        return [list(self.ppc_database['chemical_name_'+language])]

    def calculate_mean_dw_concentration(self, 
                                        tolerance = TOLERANCE_DEFAULT,
                                        relaxation_factor = RELAXATION_FACTOR_DEFAULT,
                                        max_iterations = MAX_ITERATIONS_DEFAULT):
        '''
        Calculates the mean concentration in drinking water for a 24 hour period
        given a groundwater concentration. 
        
        Parameters
        ----------
        tolerance: float 
            The allowable difference between the calculated and actual drinking water concentration, [-]
        relaxation_factor: float
            Used to iterate and calculate the new drinking water concentration, recommended 0.3-0.7 [-]
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme
        
        Returns
        -------
        mean_concentration_pipe_drinking_water: float
            Calculates the mean concentration in drinking water for a 24 hour period
            given a groundwater concentration.

        '''
        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor
        
        # Check if the flow rate, conditions have been set and parameters 
        # validated, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. \
            To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. \
            To set validate use .validate_input_parameters() ')

        else: 
            concentration_drinking_water = 0 #initial guess for drinking water 
            counter = 0

            while True:    

                sum_mass_segment = 0

                for segment in self.segment_list:

                    segment._calculate_mean_dw_mass_per_segment(concentration_groundwater=self.concentration_groundwater,
                                                                concentration_drinking_water = concentration_drinking_water,
                                                                _conditions_set = self._conditions_set,
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
                
            self.mean_concentration_pipe_drinking_water = concentration_pipe_drinking_water

        return concentration_pipe_drinking_water 

    def calculate_peak_dw_concentration(self, 
                                        tolerance = TOLERANCE_DEFAULT,
                                        relaxation_factor = RELAXATION_FACTOR_DEFAULT,
                                        max_iterations = MAX_ITERATIONS_DEFAULT):

        '''
        Calculates the peak (maximum) concentration in drinking water for a 
        given a stagnation period given a groundwater concentration.
        Stagnation period default of 8 hours. 
        
        Parameters
        ----------
        tolerance: float 
            The allowable difference between the calculated and actual drinking water concentration, [-]
        relaxation_factor: float
            Used to iterate and calculate the new drinking water concentration, recommended 0.3-0.7 [-]
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        Returns
        -------
        peak_concentration_pipe_drinking_water: float
            Calculates the peak (maximum) concentration in drinking water for a 
            given a stagnation period given a groundwater concentration.

        '''

        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor

        if self.stagnation_time != 8 * 60 * 60: #ah_todo write test for this
            print("Warning: the stagnation factor is only valid for a stagnation \
                  time of 8 hours. Using a different stagnation time is not advised.")

        # Check if the flow rate, conditions have been set and parameters 
        # validated, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. \
            To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. \
            To set validate use .validate_input_parameters() ')

        else: 
            concentration_drinking_water = 0.0 #initial guess for drinking water
            counter = 0

            while True:    

                sum_mass_segment = 0

                for segment in self.segment_list:
                    segment._calculate_peak_dw_mass_per_segment(concentration_groundwater=self.concentration_groundwater,
                                                                concentration_drinking_water = concentration_drinking_water,
                                                                _conditions_set = self._conditions_set,
                                                                stagnation_time = self.stagnation_time, 
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
                
            self.peak_concentration_pipe_drinking_water = concentration_pipe_drinking_water

        return concentration_pipe_drinking_water 


    def calculate_mean_allowable_gw_concentration(self, #ah_todo write test
                                        concentration_drinking_water,
                                        chemical_name,
                                        temperature_groundwater,
                                        language = 'NL',
                                        tolerance = TOLERANCE_DEFAULT,
                                        relaxation_factor = RELAXATION_FACTOR_DEFAULT,
                                        max_iterations = MAX_ITERATIONS_DEFAULT, 
                                        ):
        '''
        Calculates the mean 24 hour concentration in groundwater which would not 
        result in a drinking water concentration exceeding the drinking water
        norm. If the distribution coefficient it unknown for 
        a given chemical, no soil concentration is calculated. 
        
        #ah_todo add soil

        Parameters
        ----------
        concentration_drinking_water: float
            Concentration in the drinking water for which to calculate the mean 
            allowable groundwater concentration, g/m^3
        temperature_groundwater #ah_todo finish
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible
        tolerance: float 
            The allowable difference between the calculated and actual drinking water concentration, [-]
        relaxation_factor: float
            Used to iterate and calculate the new drinking water concentration, recommended 0.3-0.7 [-]
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        Returns
        -------
        concentration_mean_allowable_groundwater: float
            Mean concentration in groundwater which would would not result in 
            a mean daily (24 hours) concentration in drinking water exceeding 
            the drinking water norm, g/m3.
                    
        '''
        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor

        # Check if the flow rate, conditions have been set and parameters 
        # validated, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. \
            To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. \
            To set validate use .validate_input_parameters() ')
        
        else: 
            self._fetch_chemical_database(chemical_name=chemical_name, 
                                          suppress_print = True,
                                           language=language )

            # calculate initial guess for gw concentration
            sum_KDA_d = 0
            for segment in self.segment_list:
                # calculate the sum of the Kpw * DP * SA / d for all pipe segments
                log_Dp_ref = segment._calculate_ref_logD(chemical_group_number=self.chemical_group_number,
                            molecular_weight=self.molecular_weight)
                log_Kpw_ref = segment._calculate_ref_logK(chemical_group_number=self.chemical_group_number,
                            log_octanol_water_partitioning_coefficient=self.log_octanol_water_partitioning_coefficient)
                
                sum_KDA_d_segment = ( 10 ** log_Dp_ref * 10 ** log_Kpw_ref * segment.permeation_surface_area 
                                    / segment.diffusion_path_length )

                sum_KDA_d += sum_KDA_d_segment

            concentration_groundwater = (concentration_drinking_water * (1
                                         + self.flow_rate * segment.ASSESSMENT_FACTOR_GROUNDWATER ) 
                                            / sum_KDA_d )
            counter = 0

            while True:
                self.set_conditions(chemical_name=chemical_name, 
                                            temperature_groundwater=temperature_groundwater, 
                                            concentration_groundwater=concentration_groundwater,
                                            flow_rate = self.flow_rate,
                                            suppress_print = True, )
                               
                sum_mass_segment = 0

                # mass of chemical in pipe water to meet drinking water norm
                mass_drinkingwater_norm = concentration_drinking_water * self.flow_rate

                for segment in self.segment_list:
                    segment._calculate_mean_dw_mass_per_segment(concentration_groundwater=self.concentration_groundwater,
                                                                concentration_drinking_water = concentration_drinking_water,
                                                                    _conditions_set = self._conditions_set,
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
        
        self.concentration_mean_allowable_groundwater = concentration_groundwater

        return concentration_groundwater 


    def calculate_peak_allowable_gw_concentration(self, #ah_todo write test
                                    concentration_drinking_water,
                                    chemical_name,
                                    temperature_groundwater,
                                    language = 'NL',
                                    stagnation_time = 8 * 60 * 60,
                                    tolerance = TOLERANCE_DEFAULT,
                                    relaxation_factor = RELAXATION_FACTOR_DEFAULT,
                                    max_iterations = MAX_ITERATIONS_DEFAULT

                                    ):
        '''
        Calculates the peak (maximum) concentration in groundwater water for a 
        given a stagnation period that would not result in a peak concentration 
        in drinking water exceeding the drinking water norm for each pipe segment.
        Stagnation period default of 8 hours. If the distribution 
        coefficient it unknown for a given chemical, no soil concentration is 
        calculated.

        #ah_todo add soil

        Parameters
        ----------
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible

        temperature_groundwater  #ah_todo finish      
        stagnation_time: float
            Time in seconds which water in pipe is stagnant, unit of seconds. The 
            stagnation factor is only valid for a stagnation time of 8 hours 
            (28800 seconds), therefore using another other stagnation time is not advised.
        tolerance: float 
            The allowable difference between the calculated and actual drinking 
            water concentration, [-]
        relaxation_factor: float
            Used to iterate and calculate the new drinking water concentration, recommended 0.3-0.7 [-]
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        Returns
        -------
        concentration_peak_allowable_groundwater: float
            Concentration in groundwater which, after a stagnation period, 
            would not result in a peak concentration in drinking water exceeding 
            the drinking water norm, g/m3.
            
        '''
        self.max_iterations = int(max_iterations)
        self.tolerance = tolerance
        self.relaxation_factor = relaxation_factor
        self.stagnation_time = stagnation_time

        if stagnation_time != 8 * 60 * 60:
            print("Warning: the stagnation factor is only valid for a stagnation \
                  time of 8 hours. Using a different stagnation time is not advised.")

        # Check if the flow rate, conditions have been set and parameters 
        # validated, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. \
            To set flow rate use .set_flow_rate()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. \
            To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. \
            To set validate use .validate_input_parameters() ')

        else: 

            self._fetch_chemical_database(chemical_name=chemical_name, 
                                          suppress_print = True, 
                                          language=language)

            # calculate initial guess for gw concentration
            sum_KDA_d = 0
            for segment in self.segment_list:
                # calculate the sum of the Kpw * DP * SA *f_stag / d for all pipe segments
                log_Dp_ref = segment._calculate_ref_logD(chemical_group_number=self.chemical_group_number,
                            molecular_weight=self.molecular_weight)
                log_Kpw_ref = segment._calculate_ref_logK(chemical_group_number=self.chemical_group_number,
                            log_octanol_water_partitioning_coefficient=self.log_octanol_water_partitioning_coefficient)
                
                #stagnation factor with reference values for LogDp and LogKpw
                stagnation_factor = 10 ** max((((log_Dp_ref + 12.5) / 2 + 
                                    log_Kpw_ref) * 0.73611 + 
                                    -1.03574 ), 0)            

                sum_KDA_d_segment = ( 10 ** log_Dp_ref * 10 ** log_Kpw_ref * segment.permeation_surface_area 
                                    * stagnation_factor
                                    / segment.diffusion_path_length )

                sum_KDA_d += sum_KDA_d_segment

            # initial guess concentration in groundwater
            concentration_groundwater = concentration_drinking_water * (1 
                                        + self.total_volume * segment.ASSESSMENT_FACTOR_GROUNDWATER 
                                        / self.stagnation_time / sum_KDA_d)
           
            counter = 0

            while True:
                self.set_conditions(chemical_name=chemical_name, 
                                            temperature_groundwater=temperature_groundwater, 
                                            concentration_groundwater=concentration_groundwater, 
                                            flow_rate = self.flow_rate,
                                            suppress_print = True,)               
                sum_mass_segment = 0

                # mass of chemical in pipe water to meet drinking water norm
                mass_drinkingwater_norm = (concentration_drinking_water * self.total_volume)
                
                for segment in self.segment_list:
                    segment._calculate_peak_dw_mass_per_segment(concentration_groundwater=self.concentration_groundwater,
                                                                concentration_drinking_water = concentration_drinking_water,
                                                                    _conditions_set = self._conditions_set,
                                                                    stagnation_time = stagnation_time,
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

        self.concentration_peak_allowable_groundwater = concentration_groundwater

        return concentration_groundwater 


