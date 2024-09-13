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

from pipepermcalc.segment import *
from pipepermcalc.database import DATA_DIR 
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")

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
    _concentration_groundwater_set: Boolean
        Default False, True when the flow rate has been set.
    _is_validated: Boolean
        Default False, True when the flow rate has been set.
    total_volume: float
        Total volume of the pipe, summed from the pipe segments, m3
    total_length': float
        Total length of the pipe, summed from the pipe segments, m
    flow_rate: float
        Flow rate through pipe, m3/day.
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
    _Kd_known: Boolean
        True when a distribution coefficient (Kd) for the chemical is known. 
        Kd is needed to convert from soil to groundwater concentration.
    concentration_groundwater: float
        Concentration of the given chemical in groundwater, g/m3.
    tolerance: float 
        The degree of acceptable error in the accuracy of the calculation, 
        calculated as the difference between the calculated and actual drinking 
        water concentration, default value of 0.01 (1%), [-].
    max_iterations: int
        Maximum number of iterations allowed in the optimization scheme. 
    temperature_groundwater: float
        Temperature of the groundwater, degrees Celcius.
    stagnation_time: float
        Time in seconds which water in pipe is stagnant, unit of seconds. The 
        stagnation factor is only valid for a stagnation time of 8 hours 
        (28800 seconds), therefore using another other stagnation time is not advised.
    concentration_soil: float
        Concentration of the given chemical in soil, mg/kg.
    scale_factor_upper_limit: float
        Scale factor used to set the upper limit of the bounds for calculating 
        the mean concentration of drinking water or groundwater. Upper limit taken as the 
        concentration of groundwater (solving for drinking water concentration) 
        or solubility (solving for groundwater concentration) multiplied by the 
        scale factor. Default value of 0.999 
    scale_factor_lower_limit: float    
        Scale factor used to set the upper limit of the bounds for calculating 
        the mean concentration of drinking water or groundwater. Lower limit taken as the 
        concentration of groundwater (solving for drinking water concentration) 
        or solubility (solving for groundwater concentration) multiplied by the 
        scale factor. Default value of 0.0001.  
    ASSESSMENT_FACTOR_GROUNDWATER: float 
        Factor used to correct calculations for observations in actual pipe 
        permeation. Permeation of PE house connections in groundwater = 3, 
        other pipe materials = 1. See section 7.2 in KWR 2016.056
    ASSESSMENT_FACTOR_SOIL: float
        Factor used to correct calculations for observations in actual pipe 
        permeation. All pipe materials = 1.
            
    Note
    ----
    All parameters are in SI units: m, m2, g/m3 (equivalent to mg/L), seconds. 
    Except for the flow rate, which is given in m3/day.

    '''

    #Constants for iterative calculations,
    TOLERANCE_DEFAULT = 0.01
    SCALE_FACTOR_UPPER_LIMIT = 0.999
    SCALE_FACTOR_LOWER_LIMIT = 0.0001
    MAX_ITERATIONS_DEFAULT = 1000
    TEMPERATURE_GROUNDWATER_DEFAULT = 12 # degrees C
    STAGNATION_TIME_DEFAULT = 8 * 60 * 60 # 8 hours in seconds
    ASSESSMENT_FACTOR_GROUNDWATER = 3
    ASSESSMENT_FACTOR_SOIL = 1
    # FLOW_RATE_DEFAULT = 0.5/ 24 / 60 / 60 # 0.5 m3/day in seconds 
    # @martin, do we want a default flow rate?

    ppc_database = pd.read_csv(DATA_DIR + '/ppc_database.csv',  skiprows=[1, 2] ) 
    
    # drop chemicals for which the partitiong and/or diffusion ccoefficient cannot be calculated
    ppc_database = ppc_database.dropna(subset=['molecular_weight', 'solubility'])

    #dictionary to check input parameter value and dtype
    parameter_validation_dictionary = \
        {
        'name':     {'value_dtype': [str]}, 
        'material': {'str_options': ['PE40', 'PE80', 'SBR', 'EPDM', 'PVC'],
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
        'stagnation_time': {'min_value': 0, 
                            'value_dtype': [float, int]},  
        'flow_rate': {'min_value': 0, 
                    'value_dtype': [float, int]},  
        'concentration_soil': {'min_value': 0, 
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
        Attributes of the class added, default values of False added for the 
        different conditions (flow rate, concentration groundwater, validation etc).

        Attributes
        ----------
        segment_list: list
            list of the segments objects
        _conditions_set: Boolean
            Default False, True when the groundwater conditions have been set.
        _flow_rate_set: Boolean
            Default False, True when the flow rate has been set.
        _concentration_groundwater_set: Boolean
            Default False, True when the flow rate has been set.
        _is_validated: Boolean
            Default False, True when the flow rate has been set.
        total_volume: float
            Total volume of the pipe, summed from the pipe segments, m3
        total_length': float
            Total length of the pipe, summed from the pipe segments, m
                
        ''' 

        self.segment_list = segment_list
        self._conditions_set = False
        self._flow_rate_set = False
        self._concentration_groundwater_set = False
        self._is_validated = False
        self._set_permeation_direction = False
        self._check_not_all_PVC = False

        sum_total_volume = 0
        sum_total_length = 0
        for segment in segment_list:
            sum_total_length += segment.length
            sum_total_volume += segment.volume

        self.total_length = sum_total_length
        self.total_volume = sum_total_volume
    

    def _validate_object(self, 
                         check_object):
        ''' Check that the input parameters are valid values and types for the 
        input object, uses parameter_validation_dictionary to verify parameters.
        
        Parameters
        ----------
        check_object: object
            Segment or Pipe object for which attributes must be validated.
        '''

        for k, v in self.parameter_validation_dictionary.items():
            
            if hasattr(check_object, k):
                if (k == 'flow_rate' or k=='concentration_groundwater' or k=='concentration_soil') and getattr(check_object, k) is None:
                    pass
                else:
                    if type(getattr(check_object, k)) not in v['value_dtype']:
                        raise ValueError(f"Invalid value ~{getattr(check_object, k)}~ for parameter {k}. Input value should be a {v['value_dtype']}.")
                    if 'min_value' in v.keys():
                        if getattr(check_object, k) <= v['min_value']:
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

        #check if conditions are set
        if self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. To set pipe conditions use .set_conditions() ')
        else: 
            #validate the segment attributes
            for segment in self.segment_list:
                self._validate_object(segment)
                if segment.permeation_direction == 'perpendicular':
                    self._set_permeation_direction = True 
        
            #validate the pipe attributes
            self._validate_object(self)

            if self._set_permeation_direction is False:
                raise ValueError('Error, a pipe must consist of a least one segment with permeation perpendicular to the flow direction. Segments with permeation parallel to the flow direction are meant to join pipe segments (e.g. o-rings) with perpendicular permeation.')
            else:
                self._is_validated=True


    def _fuzzy_min_score(self, 
                         chemical_name): 
        """
        This function calculates the minimum score required for a valid
        match in fuzzywuzzy's extractOne function. The minimum score depends
        on the length of 's' and is calculated based on the string lengths and
        scores in the DEFAULT_MINSCORES dictionary. 
        From Vincent Post. 

        Parameters
        ----------
        chemical_name : str
            String for which the minimum score must be determined, in our case 
            the chemical name

        Returns
        -------
        result : float
            The minimum score for 'chemical_name'.
        """
        DEFAULT_FUZZY_MINSCORES = {1: 100, 3: 100, 4: 90, 5: 85, 6: 80, 8: 75}

        xp = list(DEFAULT_FUZZY_MINSCORES.keys()) 
        fp = [v for v in DEFAULT_FUZZY_MINSCORES.values()]
        # Use the interp function from NumPy. By default this function
        # yields fp[0] for x < xp[0] and fp[-1] for x > xp[-1]
        return np.interp(len(chemical_name), xp, fp)


    def _extract_matching_chemical_name(self, 
                                        chemical_name, 
                                        database):
        
        ''' 
        Search and extract the highest matching chemical name from the 
        database for the given input.

        Parameters
        ----------
        chemical_name: str
            String for which the minimum score/highest matching chemical name 
            must be determined.
        database: pandas df
            Dataframe of the database of chemical information, including chemical 
            name, CAS number, molecular weight, solubility, LKow, Kd etc. 

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
        chemical choice and corresponding pipe material. 

        Parameters
        ----------
        chemical_name: str
            Name of the chemical for which to calculate the permeation
        suppress_print: Boolean
            Suppress printing the chemical name and matching name, e.g. in 
            loop calculations
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible
        '''
        
        database = list(self.ppc_database['chemical_name_'+language])
        
        matching_chemical_name = self._extract_matching_chemical_name(chemical_name=chemical_name, 
                                             database=database)
        
        if suppress_print:
            pass
        else:
            print("Input chemical name: ", str(chemical_name), "- Matched chemical name: ", str(matching_chemical_name))

        df = self.ppc_database.loc[self.ppc_database['chemical_name_'+language] == matching_chemical_name]
        pipe_permeability_dict = df.to_dict('records')[0]
        self.chemical_information = pipe_permeability_dict

        #assign dict items as attribute of class
        for k, v in pipe_permeability_dict.items():
            setattr(self, k, v)

    def _groundwater_to_soil(self):
        ''' 
        Calculate the concentration in soil (mg/kg) given the concentration in 
        groundwater (g/m3 == mg/L)
        '''
        concentration_soil = (10 ** self.log_distribution_coefficient * self.concentration_groundwater 
                              * self.ASSESSMENT_FACTOR_SOIL / self.ASSESSMENT_FACTOR_GROUNDWATER)

        return concentration_soil

    def set_conditions(self,
                    chemical_name,                                    
                    concentration_groundwater=None,
                    concentration_soil=None,
                    flow_rate=None,
                    concentration_drinking_water=None,
                    temperature_groundwater= TEMPERATURE_GROUNDWATER_DEFAULT, 
                    stagnation_time = STAGNATION_TIME_DEFAULT,
                    suppress_print = False, 
                    suppress_warning = False,
                    language = 'NL'
                    ):
        ''' 
        Specifies the chemical of interest, concentration and temperature in the 
        groundwater and returns the parameters as attributes of the class. 
        If the concentration of groundwater is given, or the soil concentration 
        and Kd are known, the diffusion and permeation parameters are calculated 
        for the pipe segment(s). 

        Parameters
        ----------
        chemical_name: string
            Name of the chemical for which to calculate the permeation
        concentration_groundwater: float
            Concentration of the given chemical in groundwater, g/m3
        concentration_soil: float
            Concentration of the given chemical in soil, mg/kg.
        flow_rate: float
            Flow rate through pipe, m3/day.
        concentration_drinking_water: float
            Concentration of given chemical in drinking water pipe, g/m3. If no 
            value given, concentration is assigned the drinking water norm value.
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
        self.temperature_groundwater = temperature_groundwater
        self.stagnation_time = stagnation_time

        if (language == 'EN') or (language == 'NL'):
            self.language = language
        else: 
            raise ValueError('Error, language must be one of "EN" or "NL" ')      

        self._conditions_set = True

        self.flow_rate = flow_rate
        if flow_rate is not None:
            self._flow_rate_set = True

        self._fetch_chemical_database(chemical_name=self.chemical_name, 
                                        suppress_print=suppress_print, 
                                        language=language)


        #check if there is a known distribution coefficient
        if np.isnan(self.log_distribution_coefficient):
            self._Kd_known = False
        else: self._Kd_known = True

        self.concentration_groundwater = concentration_groundwater
        self.concentration_soil = concentration_soil

        if (concentration_groundwater is None) and (concentration_soil is None):
            pass #values already assigned, are None
        if (concentration_groundwater is not None) and (concentration_soil is None) and (self._Kd_known): 
            self.concentration_soil = self._groundwater_to_soil()
        if (concentration_groundwater is None) and (concentration_soil is not None) and (self._Kd_known): 
            self.concentration_groundwater = ((self.concentration_soil * self.ASSESSMENT_FACTOR_GROUNDWATER) 
                                                / ( 10 ** self.log_distribution_coefficient * self.ASSESSMENT_FACTOR_SOIL ))
        if (concentration_groundwater is not None) and (concentration_soil is not None) and (self._Kd_known):
            # @martin, take the gw concentration over the given soil concentration?
            # or check the Kd? 
            print('Warning, both concentration_groundwater and concentration_soil given. Only using concentration_groundwater')
            self.concentration_soil = self._groundwater_to_soil()

        if self.concentration_groundwater is not None:
            self._concentration_groundwater_set = True

        # The default value for the concentration_drinking_water is the drinking water norm
        if concentration_drinking_water is None:
            if suppress_warning:
                pass
            else:
                print(f'Warning, no drinking water concentration was defined so the drinking water concentration has been set to the norm value: ({self.Drinking_water_norm}) g/m3.')
            self.concentration_drinking_water = self.Drinking_water_norm

        else: 
            self.concentration_drinking_water = concentration_drinking_water

        if self.concentration_groundwater is not None: 
            # self.validate_input_parameters()
            for segment in self.segment_list:          
                segment._calculate_pipe_K_D(pipe = self, 
                                            _conditions_set=self._conditions_set, )


    def view_database_chemical_names(self, 
                                     language='NL'):
        '''
        Function to view a list of the possible chemical names from the database.

        Parameters
        ----------
        language: str
            Language fo the chemical name to search for, default is Dutch ('NL'), 
            English ('EN') also possible
        '''

        return list(self.ppc_database['chemical_name_'+language])

    def _check_max_terations_tolerance(self,
                                       max_iterations, 
                                       tolerance):
        ''' 
        Check the input values for max_iterations and the tolerance.

        Parameters
        ----------
        tolerance: float 
            The degree of acceptable error in the accuracy of the calculation, 
            calculated as the difference between the calculated and actual drinking 
            water concentration, default value of 0.01 (1%), [-].
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme
        '''

        if max_iterations <= 0:
            raise ValueError('Error, max_iterations must be > 0')
        else: 
            self.max_iterations = int(max_iterations)

        if (tolerance <= 0) or (tolerance > 1):
            raise ValueError('Error, tolerance must be between 0 and 1')
        else: 
            self.tolerance = tolerance

    def calculate_mean_dw_concentration(self, 
                                        tolerance = TOLERANCE_DEFAULT,
                                        max_iterations = MAX_ITERATIONS_DEFAULT,
                                        debug=False,):
        '''
        Calculates the mean concentration in drinking water for a 24 hour period
        given a groundwater concentration. 
        
        Parameters
        ----------
        tolerance: float 
            The degree of acceptable error in the accuracy of the calculation, 
            calculated as the difference between the calculated and actual drinking 
            water concentration, default value of 0.01 (1%), [-].
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        Returns
        -------
        mean_concentration_pipe_drinking_water: float
            Calculates the mean concentration in drinking water for a 24 hour period
            given a groundwater concentration.

        '''
        self._check_max_terations_tolerance(max_iterations=max_iterations, 
                                       tolerance=tolerance)
        
        # Check if the flow rate, conditions have been set and parameters 
        # validated, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. Input the flow rate in .set_conditions()')

        if self._concentration_groundwater_set is False: 
            raise ValueError('Error, the groundwater concentration has not been set. Input the groundwater concentration in .set_conditions()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. To set validate use .validate_input_parameters() ')

        else: 
            counter = 0
            concentration_drinking_water_n_plus_1 = 0 #initial guess for drinking water 
            lower_limit = 0 # initial value for the lower limit
            upper_limit = self.concentration_groundwater # initial value for the upper limit
            goodness_fit_list = [0] # initial list of goodness_fit values
            while True: 

                # counter for the number of loops in the while statement, used to prevent infinite looping
                counter +=1
                # set the drinking water concentration to be the updated one from the last loop
                concentration_drinking_water_n_min_1 = concentration_drinking_water_n_plus_1
                # initalize summed mass of chemical in segments to 0
                sum_mass_segment = 0
                # calculate the mean dw mass per segment, and sum all segments for the pipe
                for segment in self.segment_list:
                    segment._calculate_mean_dw_mass_per_segment(pipe=self, 
                                            concentration_drinking_water=concentration_drinking_water_n_min_1,
                                            concentration_groundwater=self.concentration_groundwater,)
                    sum_mass_segment += segment.mass_chemical_drinkwater

                #calculate the dw concentration from the summed segments
                concentration_drinking_water_n = (sum_mass_segment / 
                                                self.flow_rate ) 
                
                # goodness_fit used to test if the calculated dw concentration is within the tolerance value for a correct solution
                goodness_fit = abs(1 - concentration_drinking_water_n_min_1 
                                    / concentration_drinking_water_n)
                
                # check if goodness_fit value meets the allowed tolerance value, 
                # break loop is maximum iterations is exceeded
                # if not calculate the new dw concentration
                if goodness_fit <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                else:
                    goodness_fit_list.append(goodness_fit)
                    
                    # For the first two iterations use two set concentrations to begin the search for the solution, 
                    # after the 2nd iteration, search for solution by half interval search
                    if counter == 1: 
                        concentration_drinking_water_n_plus_1 = self.concentration_groundwater * self.SCALE_FACTOR_UPPER_LIMIT 
                    if counter == 2:
                        concentration_drinking_water_n_plus_1 = self.concentration_groundwater * self.SCALE_FACTOR_LOWER_LIMIT
                    if counter >2:
                        if (goodness_fit < goodness_fit_list[counter-1]) or (concentration_drinking_water_n > self.concentration_groundwater):
                            lower_limit = concentration_drinking_water_n_min_1
                            concentration_drinking_water_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                        else:
                            upper_limit = concentration_drinking_water_n_min_1
                            concentration_drinking_water_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                    if debug: 
                        print(concentration_drinking_water_n_min_1, concentration_drinking_water_n_plus_1, goodness_fit, lower_limit, upper_limit)

            # assign the drinking water concentration to be the concentration calculated in the loop                            
            self.concentration_drinking_water = concentration_drinking_water_n
            if concentration_drinking_water_n > self.solubility:
                print(f'Warning, the calculated drinking water concentration ({concentration_drinking_water_n}) is above the solubility limit, {self.solubility}.')

        return concentration_drinking_water_n 


    def calculate_peak_dw_concentration(self, 
                                        tolerance = TOLERANCE_DEFAULT,
                                        max_iterations = MAX_ITERATIONS_DEFAULT,
                                        debug=False):

        '''
        Calculates the peak (maximum) concentration in drinking water for a 
        given a stagnation period given a groundwater concentration.
        Stagnation period default of 8 hours. 
        
        Parameters
        ----------
        tolerance: float 
            The degree of acceptable error in the accuracy of the calculation, 
            calculated as the difference between the calculated and actual drinking 
            water concentration, default value of 0.01 (1%), [-].
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme

        Returns
        -------
        peak_concentration_pipe_drinking_water: float
            Calculates the peak (maximum) concentration in drinking water for a 
            given a stagnation period given a groundwater concentration.

        '''
        self._check_max_terations_tolerance(max_iterations=max_iterations, 
                                       tolerance=tolerance, )

        if self.stagnation_time != self.STAGNATION_TIME_DEFAULT: 
            print("Warning: the stagnation factor is only valid for a stagnation time of 8 hours. Using a different stagnation time is not advised.")

        # Check if the conditions have been set and parameters 
        # validated, if not raise error

        if self._concentration_groundwater_set is False: 
            raise ValueError('Error, the groundwater concentration has not been set. Input the groundwater concentration in .set_conditions()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. To set validate use .validate_input_parameters() ')

        else: 
            counter = 0
            concentration_drinking_water_n_plus_1 = 0
            lower_limit = 0 # initial value for the lower limit
            upper_limit = self.concentration_groundwater # initial value for the upper limit
            goodness_fit_list = [0] # initial list of goodness_fit values

            while True:                
                # counter for the number of loops in the while statement, used to prevent infinite looping
                counter +=1
                # set the drinking water concentration to be the updated one from the last loop
                concentration_drinking_water_n_min_1 = concentration_drinking_water_n_plus_1
                # initalize summed mass of chemical in segments to 0
                sum_mass_segment = 0
                # calculate the mean dw mass per segment, and sum all segments for the pipe
                for segment in self.segment_list:
                    segment._calculate_peak_dw_mass_per_segment(pipe=self, 
                                            concentration_drinking_water=concentration_drinking_water_n_min_1,
                                            concentration_groundwater=self.concentration_groundwater,)

                    sum_mass_segment += segment.mass_chemical_drinkwater
                
                #calculate the dw concentration from the summed segments
                concentration_drinking_water_n = (sum_mass_segment / 
                                                self.total_volume ) 
                # goodness_fit used to test if the calculated dw concentration is within the tolerance value for a correct solution
                goodness_fit = abs(1 - concentration_drinking_water_n_min_1 / concentration_drinking_water_n)
                
                # check if goodness_fit value meets the allowed tolerance value, 
                # break loop is maximum iterations is exceeded
                # if not calculate the new dw concentration
                if goodness_fit <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                else:
                    goodness_fit_list.append(goodness_fit)
                    
                    # For the first two iterations use two set concentrations to begin the search for the solution, 
                    # after the 2nd iteration, search for solution by half interval search
                    if counter == 1:
                        concentration_drinking_water_n_plus_1 = self.concentration_groundwater * self.SCALE_FACTOR_UPPER_LIMIT
                    if counter == 2:
                        concentration_drinking_water_n_plus_1 = self.concentration_groundwater * self.SCALE_FACTOR_LOWER_LIMIT
                    if counter >2:
                        if (goodness_fit < goodness_fit_list[counter-1]) or (concentration_drinking_water_n > self.concentration_groundwater):
                            lower_limit = concentration_drinking_water_n_min_1
                            concentration_drinking_water_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                        else:
                            upper_limit = concentration_drinking_water_n_min_1
                            concentration_drinking_water_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                    if debug: 
                        print(concentration_drinking_water_n_min_1, concentration_drinking_water_n_plus_1, goodness_fit, lower_limit, upper_limit)

            # assign the drinking water concentration to be the concentration calculated in the loop                            
            self.concentration_drinking_water = concentration_drinking_water_n
            if concentration_drinking_water_n > self.solubility:
                print(f'Warning, the calculated drinking water concentration ({concentration_drinking_water_n}) is above the solubility limit, {self.solubility}.')

        return concentration_drinking_water_n


    def calculate_mean_allowable_gw_concentration(self, 
                                        tolerance = TOLERANCE_DEFAULT,
                                        max_iterations = MAX_ITERATIONS_DEFAULT, 
                                        debug=False,):
        '''
        Calculates the mean 24 hour concentration in groundwater which would not 
        result in a drinking water concentration exceeding the drinking water
        norm. If the distribution coefficient is known the soil concentration is
        also calculated. 

        Parameters
        ----------
        tolerance: float 
            The degree of acceptable error in the accuracy of the calculation, 
            calculated as the difference between the calculated and actual drinking 
            water concentration, default value of 0.01 (1%), [-].
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme
        debug: Boolean
            If True, return the groundwater concentration, goodness_fit and lower and 
            upper limits every iteration.
 
        Returns
        -------
        concentration_mean_allowable_groundwater: float
            Mean concentration in groundwater which would would not result in 
            a mean daily (24 hours) concentration in drinking water exceeding 
            the drinking water norm, g/m3.
                    
        '''
        self._check_max_terations_tolerance(max_iterations=max_iterations, 
                                       tolerance=tolerance)

        # Check if the flow rate, conditions have been set and parameters 
        # validated, if not raise error
        if self._flow_rate_set is False: 
            raise ValueError('Error, the flow rate in the pipe has not been set. To set flow rate use .set_flow_rate()')

        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. To set validate use .validate_input_parameters() ')
        
        if self.concentration_drinking_water is None:
            raise ValueError('Error, no default drinking water norm, please input a drinking water concentration using .set_conditions()')

        if self.concentration_drinking_water > self.solubility:
            raise ValueError('Error, the drinking water concentration given or the default drinking water norm is higher than the solubility of the chemical. Input a lower drinking water concentration using .set_conditions()')

        else: 
            self._fetch_chemical_database(chemical_name=self.chemical_name, 
                                          suppress_print = True,
                                           language=self.language )

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

                # initial guess concentration in groundwater
                concentration_groundwater_n_plus_1 = (self.concentration_drinking_water * (1
                                         + self.flow_rate * self.ASSESSMENT_FACTOR_GROUNDWATER ) 
                                            / sum_KDA_d ) * 24 * 60 * 60 
            
            counter = 0
            lower_limit = self.concentration_drinking_water # initial value for the lower limit
            upper_limit = self.solubility # initial value for the upper limit
            goodness_fit_list = [0] # initial list of goodness_fit values

            while True:
                # counter for the number of loops in the while statement, used to prevent infinite looping
                counter +=1                
                # set the groundwater concentration to be the updated one from the last loop
                concentration_groundwater_n_min_1 = concentration_groundwater_n_plus_1

                # set the conditions for the pipe with the updated groundwater concentration
                self.set_conditions(chemical_name=self.chemical_name,                                    
                    concentration_groundwater=concentration_groundwater_n_min_1,
                    flow_rate=self.flow_rate,
                    concentration_drinking_water=self.concentration_drinking_water, 
                    temperature_groundwater=self.temperature_groundwater, 
                    stagnation_time = self.stagnation_time,
                    suppress_print = True, 
                    language = self.language)

                # initalize summed mass of chemical in segments to 0
                sum_mass_segment = 0

                # mass of chemical in pipe water to meet drinking water norm
                mass_drinkingwater_norm = (self.concentration_drinking_water * self.flow_rate)

                # calculate the mean dw mass per segment, and sum all segments for the pipe
                for segment in self.segment_list:
                    segment._calculate_mean_dw_mass_per_segment(pipe=self, 
                                            concentration_drinking_water=self.concentration_drinking_water,
                                            concentration_groundwater=self.concentration_groundwater,)

                    sum_mass_segment += segment.mass_chemical_drinkwater
                
                # goodness_fit used to test if the calculated dw concentration is within the tolerance value for a correct solution
                goodness_fit = abs(1 - mass_drinkingwater_norm / sum_mass_segment)

                # check if goodness_fit value meets the allowed tolerance value, 
                # break loop is maximum iterations is exceeded
                # if not calculate the new gw concentration
                if goodness_fit <= tolerance:
                    if debug:
                        print(concentration_groundwater_n_min_1, concentration_groundwater_n_plus_1, goodness_fit, lower_limit, upper_limit) #for debugging
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                elif upper_limit == lower_limit:
                    print('No solution found, lower_limit = upper_limit. Groundwater concentration necesary to satisfy the given drinking water concentration may be above the solubility limit.')
                    break
                else:
                    goodness_fit_list.append(goodness_fit)

                    # For the first two iterations use two set concentrations to 
                    # begin the search for the solution, after the 2nd iteration
                    # search for solution by half interval search
                    if counter == 1:
                        concentration_groundwater_n_plus_1 = self.solubility * self.SCALE_FACTOR_UPPER_LIMIT
                    if counter == 2:
                        concentration_groundwater_n_plus_1 = 0
                    if counter >2:
                        if (goodness_fit < goodness_fit_list[counter-1]) or (concentration_groundwater_n_plus_1 < self.concentration_drinking_water):
                            lower_limit = concentration_groundwater_n_min_1
                            concentration_groundwater_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                        else:
                            upper_limit = concentration_groundwater_n_min_1
                            concentration_groundwater_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                    if debug: 
                        print(concentration_groundwater_n_min_1, concentration_groundwater_n_plus_1, goodness_fit, lower_limit, upper_limit) #for debugging

        # assign the groundwater concentration to be the concentration calculated in the loop                            
        self.concentration_groundwater = concentration_groundwater_n_min_1

        if concentration_groundwater_n_min_1 > self.solubility:
            print(f'Warning, the calculated drinking water concentration ({concentration_groundwater_n_min_1}) is above the solubility limit, {self.solubility}.')

        if self._Kd_known: 
            self.concentration_soil = self._groundwater_to_soil()
        else: 
            self.concentration_soil = 'No known distribution coefficient to calculate soil concentration'

        return concentration_groundwater_n_min_1 


    def calculate_peak_allowable_gw_concentration(self, 
                                    tolerance = TOLERANCE_DEFAULT,
                                    max_iterations = MAX_ITERATIONS_DEFAULT,
                                    debug=False,):
        '''
        Calculates the peak (maximum) concentration in groundwater water for a 
        given a stagnation period that would not result in a peak concentration 
        in drinking water exceeding the drinking water norm for each pipe segment.
        Stagnation period default of 8 hours. If the distribution coefficient is 
        known the soil concentration is also calculated. 

        Parameters
        ----------
        tolerance: float 
            The degree of acceptable error in the accuracy of the calculation, 
            calculated as the difference between the calculated and actual drinking 
            water concentration, default value of 0.01 (1%), [-].
        max_iterations: int
            Maximum number of iterations allowed in the optimization scheme
        debug: Boolean
            If True, return the groundwater concentration, goodness_fit and lower and 
            upper limits every iteration.

            
        Returns
        -------
        concentration_peak_allowable_groundwater: float
            Concentration in groundwater which, after a stagnation period, 
            would not result in a peak concentration in drinking water exceeding 
            the drinking water norm, g/m3.
            
        '''
        self._check_max_terations_tolerance(max_iterations=max_iterations, 
                                       tolerance=tolerance)

        if self.stagnation_time != self.STAGNATION_TIME_DEFAULT:
            print("Warning: the stagnation factor is only valid for a stagnation time of 8 hours. Using a different stagnation time is not advised.")

        # Check if the conditions have been set and parameters validated, if not raise error
        elif self._conditions_set is False:
            raise ValueError('Error, the pipe conditions must first be set. To set pipe conditions use .set_conditions() ')

        elif self._is_validated is False: 
            raise ValueError('Error, the input parameters must first be validated. To set validate use .validate_input_parameters() ')

        if self.concentration_drinking_water > self.solubility:
            raise ValueError('Error, the drinking water concentration given or the default drinking water norm is higher than the solubility of the chemical. Input a lower drinking water concentration using .set_conditions()')

        else: 

            self._fetch_chemical_database(chemical_name=self.chemical_name, 
                                          suppress_print = True, 
                                          language=self.language)

            # calculate initial guess for gw concentration
            sum_KDA_d = 0
            for segment in self.segment_list:
                # calculate the sum of the Kpw * DP * SA *f_stag / d for all pipe segments
                log_Dp_ref = segment._calculate_ref_logD(chemical_group_number=self.chemical_group_number,
                            molecular_weight=self.molecular_weight)
                log_Kpw_ref = segment._calculate_ref_logK(chemical_group_number=self.chemical_group_number,
                            log_octanol_water_partitioning_coefficient=self.log_octanol_water_partitioning_coefficient)
                
                #stagnation factor with reference values for LogDp and LogKpw
                # See equation 6-10 in KWR 2016.056. Note equation description in 
                # report is given in different format. 
                stagnation_factor = 10 ** max((((log_Dp_ref + 12.5) / 2 + 
                                    log_Kpw_ref) * 0.73611 + 
                                    -1.03574 ), 0)            

                sum_KDA_d_segment = ( 10 ** log_Dp_ref * 10 ** log_Kpw_ref * segment.permeation_surface_area 
                                    * stagnation_factor / segment.diffusion_path_length )

                sum_KDA_d += sum_KDA_d_segment

            # initial guess concentration in groundwater
            concentration_groundwater_n_plus_1 = self.concentration_drinking_water * (1 
                                        + self.total_volume * self.ASSESSMENT_FACTOR_GROUNDWATER 
                                        / self.stagnation_time / sum_KDA_d) 
            
            counter = 0
            lower_limit = self.concentration_drinking_water # initial value for the lower limit
            upper_limit = self.solubility # initial value for the upper limit
            goodness_fit_list = [0] # initial list of goodness_fit values

            while True: 
                # counter for the number of loops in the while statement, used to prevent infinite looping
                counter +=1                
                # set the groundwater concentration to be the updated one from the last loop
                concentration_groundwater_n_min_1 = concentration_groundwater_n_plus_1

                # set the conditions for the pipe with the updated groundwater concentration
                self.set_conditions(chemical_name=self.chemical_name,                                    
                    concentration_groundwater=concentration_groundwater_n_min_1,
                    flow_rate=self.flow_rate,
                    concentration_drinking_water=self.concentration_drinking_water, 
                    temperature_groundwater=self.temperature_groundwater, 
                    stagnation_time = self.stagnation_time,
                    suppress_print = True, 
                    language = self.language)
                
                # initalize summed mass of chemical in segments to 0
                sum_mass_segment = 0

                # mass of chemical in pipe water to meet drinking water norm
                mass_drinkingwater_norm = (self.concentration_drinking_water * self.total_volume) 

                # calculate the mean dw mass per segment, and sum all segments for the pipe
                for segment in self.segment_list:
                    segment._calculate_peak_dw_mass_per_segment(pipe=self, 
                                            concentration_drinking_water=self.concentration_drinking_water,
                                            concentration_groundwater=self.concentration_groundwater,)

                    sum_mass_segment += segment.mass_chemical_drinkwater

                # goodness_fit used to test if the calculated dw concentration is within the tolerance value for a correct solution
                goodness_fit = abs(1 - mass_drinkingwater_norm / sum_mass_segment)

                # check if goodness_fit value meets the allowed tolerance value, 
                # break loop is maximum iterations is exceeded
                # if not calculate the new gw concentration
                if goodness_fit <= tolerance:
                    break
                elif counter > max_iterations:
                    print('Max iterations exceeded')
                    break
                elif upper_limit == lower_limit:
                    print('No solution found, lower_limit = upper_limit. Groundwater concentration necesary to satisfy the given drinking water concentration may be above the solubility limit.')
                    break
                else:
                    goodness_fit_list.append(goodness_fit)

                    # two initial guesses to compare the goodness of fit
                    if counter == 1:
                        concentration_groundwater_n_plus_1 = self.solubility * self.SCALE_FACTOR_UPPER_LIMIT
                    if counter == 2:
                        concentration_groundwater_n_plus_1 = 0
                    if counter >2:
                        if (goodness_fit < goodness_fit_list[counter-1]) or (concentration_groundwater_n_plus_1 < self.concentration_drinking_water):
                            lower_limit = concentration_groundwater_n_min_1
                            concentration_groundwater_n_plus_1 = lower_limit + (upper_limit -lower_limit)/2
                        else:
                            upper_limit = concentration_groundwater_n_min_1
                            concentration_groundwater_n_plus_1 = lower_limit - (upper_limit -lower_limit)/2
                    if debug: 
                        print(concentration_groundwater_n_min_1, concentration_groundwater_n_plus_1, goodness_fit, lower_limit, upper_limit) #for debugging
        
        # assign the groundwater concentration to be the concentration calculated in the loop
        self.concentration_groundwater = concentration_groundwater_n_min_1
        
        if concentration_groundwater_n_min_1 > self.solubility:
            print(f'Warning, the calculated drinking water concentration ({concentration_groundwater_n_min_1}) is above the solubility limit, {self.solubility}.')

        if self._Kd_known: 
            self.concentration_soil = self._groundwater_to_soil()
        else: 
            self.concentration_soil = 'No known distribution coefficient to calculate soil concentration'

        return concentration_groundwater_n_min_1 


