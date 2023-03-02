#%% ----------------------------------------------------------attributes------------------
# A. Hockin, January 2023
# KWR 403230-003
# Pipe permeation calculator
# With Martin vd Schans, Bram Hillebrand, Lennart Brokx
#
# ------------------------------------------------------------------------------
import numpy as np
import pandas as pd

from project_path import file_path

from pipepermcalc.pipe import * 

class Segment:
    ''' 
    Segment object class to make segments of the pipe.

    Attributes
    ----------
    #ah_todo add attributes

    PARTITIONING_A_DH: float
        Coefficient for correcting the partitioning coefficient for temperature. 
        From regression analysis, a is the slope, see table 5-6 in 
        KWR 2016.056. Constant equal to 7.92169801506708. 
    PARTITIONING_B_DH: float,
        Coefficient for correcting the partitioning coefficient for temperature. 
        From regression analysis, b is the intercept, see table 5-6 in 
        KWR 2016.056. Constant equal to -17.1875608983359. 
    DIFFUSION_A_DH: float 
        Coefficient for correcting the diffusion coefficient for temperature. 
        From regression analysis, a is the slope, see table 5-6 in 
        KWR 2016.056. Constant equal to 61.8565740136974. 
    DIFFUSION_B_DH: float
        Coefficient for correcting the diffusion coefficient for temperature. 
        From regression analysis, b is the intercept, see table 5-6 in 
        KWR 2016.056. Constant equal to -78.9191401984509. 
    ASSESSMENT_FACTOR_GROUNDWATER: float 
        Factor used to correct calculations for observations in actual pipe 
        permeation. Permeation of PE house connections in groundwater = 3, 
        other pipe materials = 1. See section 7.2 in KWR 2016.056
    ASSESSMENT_FACTOR_SOIL: float
        Factor used to correct calculations for observations in actual pipe 
        permeation. All pipe materials = 1.
    PARTITIONING_A_C: float
        Constant used in the correction for the partitioning coefficent due to 
        the influence of temperature. See equation 5-20 in KWR 2016.056, for 
        partitioning a_c = 0.103965019849463.
    PARTITIONING_CREF_SW: float
        Reference concentration used in the correction for the partitioning 
        coefficent due to the influence of temperature. Ssee section 5.4.7 in 
        KWR 2016.056. For partitioning, Cref_SW = 1.0.
    DIFFUSION_A_C: float
        Constant used in the correction for the diffusion coefficent due to 
        the influence of temperature. See equation 5-18 in KWR 2016.056, for 
        diffusion a_c = 0.784077209735583.
    DIFFUSION_CREF_SW: float
        Reference concentration used in the correction for the diffusion 
        coefficent due to the influence of temperature. Ssee section 5.4.6 in 
        KWR 2016.056. For partitioning, Cref_SW = 0.5.

    name: string
        name of the pipe segment
    material: enum?? #ah_todo @Bram -> set choice of materials
        e.g. PE40, PE80, PVC, EPDM, rubber etc.
    length: float
        Length of the pipe segment, meters 
    inner_diameter: float
        Inner diameter of the pipe segment, meters
    wall_thickness: float
        wall_thickness of the pipe segment, meters
    permeation_direction: string #ah_todo enum?? @Bram -> limit choice of direction
        Direction of permeation through the pipe segment. Options are 
        'perpendicular' or 'parallel'. Default permeation is perpendicular 
        to the flow direction. See schematic XX in read the docs. #ah_todo how to reference schematic
        #ah_todo limit input to enum value
    diffusion_path_length: float
        In the case of permeation parallel to the flow direction, a 
        diffusion path length is required to calculate the permeation 
        through the pipe segment. For example in the case of a pipe 
        coupling rings. If no value is given, the default value is the 
        wall_thickness, unit meters.
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
    stagnation_factor: float
        Correction for the decrease in the concentratino gradient near the 
        inner wall of the pipe during stagnation (e.g. no flow at at night)
    
    Note
    ----
    All parameters are in SI units: m, m2, g/m3 (equivalent to mg/L), seconds.


    '''

    count = 0 # count of pipe segments

    def __init__(self, 
                name=None,
                material=None,
                length=None,
                inner_diameter=None,
                wall_thickness=None,
                permeation_direction='perpendicular',
                diffusion_path_length=None,
                ):
        '''
        Creates a pipe segment with the attributes of the pipe (length, 
        wall_thickness, diameter, material etc.). 
        
        Parameters
        ----------
        name: string
            name of the pipe segment
        material: string #ah_todo enum?? @Bram -> set choice of materials
            e.g. PE40, PE80, PVC, EPDM, rubber etc.
        length: float
            Length of pipe segment, meters 
        inner_diameter: float
            Inner diameter of pipe segment, meters
        wall_thickness: float
            wall_thickness of pipe segment, meters
        permeation_direction: string #ah_todo enum?? @Bram -> limit choice of direction
            Direction of permeation through the pipe segment. Options are 
            'perpendicular' or 'parallel'. Default permeation is perpendicular 
            to the flow direction. See schematic XX in read the docs.
        diffusion_path_length: float
            In the case of permeation perpendicular to the flow direction, a 
            diffusion path length is required to calculate the permeation 
            through the pipe segment. For example in the case of a pipe 
            coupling rings. If no value is given, diffusion is assumed 
            perpendicular to the flow direction and the wall_thickness is 
            used to calculate the diffusion through the pipe segment. 
            Unit meters.
            

        '''  


        #Constants for various LogK and Log D equations
        self._PARTITIONING_A_DH = 7.92169801506708 #see table 5-6 in KWR 2016.056
        self._PARTITIONING_B_DH = -17.1875608983359 #see table 5-6 in KWR 2016.056
        self._DIFFUSION_A_DH = 61.8565740136974 #see table 5-6 in KWR 2016.056
        self._DIFFUSION_B_DH = -78.9191401984509 #see table 5-6 in KWR 2016.056
        self.ASSESSMENT_FACTOR_GROUNDWATER = 3 
        self.ASSESSMENT_FACTOR_SOIL = 1
        self.PARTITIONING_A_C = 0.103965019849463 #see equation 5-20 in KWR 2016.056
        self.PARTITIONING_CREF_SW = 1.000 #see section 5.4.7 in KWR 2016.056
        self.DIFFUSION_A_C = 0.784077209735583 #see equation 5-18 in KWR 2016.056
        self.DIFFUSION_CREF_SW = 0.5 #see section 5.4.6 in KWR 2016.056

        self.name = name
        self.material = material

        self.length = float(length)
        self.inner_diameter = float(inner_diameter)
        self.wall_thickness = float(wall_thickness)
        self.permeation_direction = permeation_direction

        if diffusion_path_length is None:
            self.diffusion_path_length = self.wall_thickness
        else:
            self.diffusion_path_length = float(diffusion_path_length)    

        outer_diameter = inner_diameter + wall_thickness
        volume = np.pi * (inner_diameter / 2) ** 2 * length
        permeation_surface_area =(np.pi * inner_diameter * length)

        if permeation_direction == 'parallel':
            volume = 0 
            permeation_surface_area = ((np.pi * ((inner_diameter + wall_thickness) ** 2 
                                                - inner_diameter ** 2) ) / 4)
                   
        self.volume = volume
        self.permeation_surface_area = permeation_surface_area
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter

    # @ah_todo revert back to csv? seperate file? 
    # From Bram, @MartinK-> suggest to implement the "named tuple" method, leave for now do at the end
    # SBR, EPDM refer to memo, ask m. meerkerk for memo #, include project number for the memo
    reference_pipe_material_dict = \
        {
        "PE40": {
            "REF_LOG_D_A": {
                1: -0.011,	
                2: -0.00629,
                3: -0.006,
                },
            "REF_LOG_D_B": {
                1: -10.688,
                2: -11.000,
                3: -11.000
                },
            "REF_LOG_K_A": {
                1: 1.097,
                2: 1.059,
                3: 0.979
                },
            "REF_LOG_K_B": {
                1: -0.689,
                2: -0.67,
                3: -1.796,
                },
        },
        "PE80": {
            "REF_LOG_D_A": {
                1: -0.011,	
                2: -0.00629,
                3: -0.00629,
                },
            "REF_LOG_D_B": {
                1: -11.188,
                2: -11.188,
                3: -11.500,
                },
            "REF_LOG_K_A": {
                1: 1.185,	
                2: 1.185,
                3: 1.231,
                },
            "REF_LOG_K_B": {
                1: -1.437,	
                2: -1.437,
                3: -2.606,
                },
        },
        "SBR": {
            "REF_LOG_D_A": {
                1: -0.011 * 0.950647410867427,	#ah_todo change these to the see notes, replace formulas in _calculate_ref_logD
                2: -0.00629 * 0.950647410867427,
                3: -0.006 * 0.950647410867427,
                },
            "REF_LOG_D_B": {
                1: -10.688 * 0.950647410867427,
                2: -11.000 * 0.950647410867427,
                3: -11.000 * 0.950647410867427,
                },
            "REF_LOG_K_A": {
                1: 1.0452,	
                2: 1.0452,
                3: 1.0452,
                },
            "REF_LOG_K_B": {
                1: -0.3686,	
                2: -0.3686,
                3: -0.3686,
                },
        },  
        "EPDM": {
            "REF_LOG_D_A": {
                1: -0.011* 0.920996123470591,	
                2: -0.00629 * 0.920996123470591,
                3: -0.006 * 0.920996123470591,
                },
            "REF_LOG_D_B": {
                1: -10.688 * 0.920996123470591,
                2: -11.000 * 0.920996123470591,
                3: -11.000 * 0.920996123470591,
                },
            "REF_LOG_K_A": {
                1: 1.0675,	
                2: 1.0675,
                3: 1.0675,
                },
            "REF_LOG_K_B": {
                1: -0.3002,	
                2: -0.3002,
                3: -0.3002,
                },
        },      
        }

    def _correct_for_temperature(self, #ah_todo @Bram, when to use or not use self in functions?
                                temperature_groundwater,
                                coefficient_name=None, 
                                a_dh=None,
                                b_dh=None,
                                ):
        '''
        Temperature correction for the partitioning and diffusion coefficients, 
        
        See table 5-3 in KWR 2016.056

        Parameters
        ----------
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

        R = 0.008314 #universal gas constant
        reference_temperature = 25 # deg. C
        dh = a_dh * np.log10(coefficient_name) + b_dh
        f_temp = dh / (R * np.log(10)) * (1 / (reference_temperature + 273) - 1 / (temperature_groundwater + 273))
        return f_temp


    def _concentration_correction(self,
                        solubility, 
                        concentration_groundwater,
                        a_c=None,
                        Cref_Sw=None):
        '''
        Correction factor for the influence of concentration on the 
        partitioning or diffusion coefficient 

        See table 5-3, equations 5-17 and 5-19 in KWR 2016.056

        Parameters
        ----------
        #ah_todo finish

        Returns
        -------
        f_conc: string
            Concentration correction factor for the partitioning or diffusion 
            coefficient
        '''

        Cg_Sw = min(concentration_groundwater / solubility, 1)
        f_conc = a_c * (Cg_Sw - Cref_Sw)

        return f_conc


    def _correct_for_age(self,):
        '''
        Age correction, none implemented yet'''

        f_age = 0.000
        return f_age


    def _calculate_ref_logK(self,
                            chemical_group_number,
                            log_octanol_water_partitioning_coefficient):
        '''Calculate the reference log K'''

        a_ref = self.reference_pipe_material_dict[self.material]['REF_LOG_K_A'][chemical_group_number]
        b_ref = self.reference_pipe_material_dict[self.material]['REF_LOG_K_B'][chemical_group_number]
        log_Kpw_ref = a_ref * log_octanol_water_partitioning_coefficient + b_ref

        return log_Kpw_ref


    def _calculate_ref_logD(self,
                            chemical_group_number,
                            molecular_weight):
        '''Calculate the reference log D based on the pipe material. A fixed 
        ratio between the log of the diffusion coefficient of PE-40 (logD_p) 
        and of SBR/EPDM (logD_s, logD_e) in m2/s is assumed for SBR and 
        EPDM for the determination of the diffusion coefficient, see memo 2022 
        "Permeatie door rubber afdichtingen van drinkwaterleidingen."  '''

        a_ref = self.reference_pipe_material_dict[self.material]['REF_LOG_D_A'][chemical_group_number]
        b_ref = self.reference_pipe_material_dict[self.material]['REF_LOG_D_B'][chemical_group_number]
        log_Dp_ref = a_ref * molecular_weight + b_ref

        return log_Dp_ref    


    def _calculate_logK(self, pipe):
        ''' 
        Calculate the LogK value for the pipe material, correct for temperature,
        concentration and age. 
        
        See table 5-3 in KWR 2016.056 for explanation of calculations
        
        '''

        # calculate reference log K plastic-water (log kpw) 
        log_Kpw_ref = self._calculate_ref_logK(chemical_group_number=pipe.chemical_group_number,
                                               log_octanol_water_partitioning_coefficient= pipe.log_octanol_water_partitioning_coefficient)

        # correct for temperature, concentration, age
        f_Ktemp = self._correct_for_temperature(coefficient_name = pipe.solubility,
                            temperature_groundwater = pipe.temperature_groundwater, 
                            a_dh = self._PARTITIONING_A_DH, 
                            b_dh = self._PARTITIONING_B_DH, )

        f_Kconc = self._concentration_correction(solubility=pipe.solubility, 
                                concentration_groundwater = pipe.concentration_groundwater,
                                a_c = self.PARTITIONING_A_C,
                                Cref_Sw = self.PARTITIONING_CREF_SW) 
        
        f_Kage = self._correct_for_age()

        # sum corrections for final Log k
        log_Kpw = log_Kpw_ref + f_Ktemp + f_Kconc + f_Kage
        
        self.log_Kpw_ref = log_Kpw_ref
        self.f_Ktemp = f_Ktemp   
        self.f_Kconc = f_Kconc
        self.log_Kpw = log_Kpw

        return log_Kpw


    def _calculate_logD(self, pipe):
        ''' 
        Calculate the LogK value for the pipe material, correct for temperature,
        concentration and age. 

        See table 5-3 in KWR 2016.056 for explanation of calculations
        
        '''
        
        # calculate reference log D plastic (log Dp) 
        log_Dp_ref = self._calculate_ref_logD(chemical_group_number=pipe.chemical_group_number,
                            molecular_weight=pipe.molecular_weight)

        # correct for temperature, concentration, age
        f_Dtemp = self._correct_for_temperature(coefficient_name =pipe.molecular_weight,
                            temperature_groundwater = pipe.temperature_groundwater, 
                            a_dh = self._DIFFUSION_A_DH, 
                            b_dh = self._DIFFUSION_B_DH,)

        f_Dconc = self._concentration_correction(
                                solubility=pipe.solubility, 
                                concentration_groundwater = pipe.concentration_groundwater,
                                a_c = self.DIFFUSION_A_C , 
                                Cref_Sw = self.DIFFUSION_CREF_SW) 
        
        f_Dage = self._correct_for_age()

        # sum corrections for final Log D
        log_Dp = log_Dp_ref + f_Dtemp + f_Dconc + f_Dage



        #ah discussed wtih Bram not to store these values, only calculate them on the fly
        self.log_Dp_ref = log_Dp_ref
        self.f_Dtemp = f_Dtemp    
        self.f_Dconc = f_Dconc
        self.log_Dp = log_Dp #m2/s

        return log_Dp

    # def _calculate_permeation_coefficient(self,):
    #     ''' Calculate the permeation coefficient for the segment'''
    #     #Permeation coefficient for plastic-water (Ppw), unit: m2/day
    #     permeation_coefficient = (24 * 60 * 60 * 
    #                                     (10 ** self.log_Dp) 
    #                                     * 10 ** self.log_Kpw)
    #     return permeation_coefficient

    def _calculate_pipe_K_D(self,
                            pipe,
                            _conditions_set, 
        ):
        '''
        Fetch the pipe and chemical information corresponding to the given pipe 
        material and chemical choice, assign as attributes of the class.

        See table 5-3 in KWR 2016.056 for explanation of calculations

        Parameters
        ----------
        pipe_material: string
            Choice of pipe material: PE40, PE80, SBR, EPDM
        '''
        # Check if the groundwater conditions have been set, if not raise error
        if _conditions_set is None:
            raise ValueError('Error, pipe conditions have not been set. \
                             To set pipe conditions use .set_conditions()')
        else:           

            # calculate log K plastic-water (log kpw) 
            self.log_Kpw = self._calculate_logK(pipe = pipe)

            # calculate log D plastic (log Dp) 
            self.log_Dp = self._calculate_logD(pipe = pipe)
            

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
                                            pipe,): 
        '''
        Calculates the mean mass in drinking water for a 24 hour period given a 
        groundwater concentration, for each pipe segment.
        
        Parameters
        ----------
        #ah_todo finish
        '''
         
        # From equation 4-7 in KWR 2016.056, but not simplifying the mass flux 
        # in equation 4-5 
        delta_c = pipe.concentration_groundwater - pipe.concentration_drinking_water

        self.mass_chemical_drinkwater = (((10 ** self.log_Dp * 10 ** self.log_Kpw)
                                          * self.permeation_surface_area 
                                          * delta_c / self.diffusion_path_length ) 
                                            / self.ASSESSMENT_FACTOR_GROUNDWATER
                                            * 24 * 60 * 60)


    def _calculate_peak_dw_mass_per_segment(self, 
                                            pipe,):
        '''
        Calculates the peak (maximum) mass in drinking water for a 
        given a stagnation period given a groundwater concentration, for each pipe segment.
        Stagnation period default of 8 hours. 
        
        Parameters
        ----------
        pipe_segment: string
            name of the pipe segment        
        stagnation_time: float
            Time in seconds which water in pipe is stagnant, unit of seconds. The 
            stagnation factor is only valid for a stagnation time of 8 hours 
            (28800 seconds), therefore using another other stagnation time is not advised.
        #ah_todo finish

        '''

        self.stagnation_factor = self._calculate_stagnation_factor()
        delta_c = pipe.concentration_groundwater - pipe.concentration_drinking_water

        # From equation 4-10 KWR 2016.056, but not simplifying the mass flux 
        # in equation 4-5 and rearranging to remove C_dw from the equation       
        self.mass_chemical_drinkwater = (((10 ** self.log_Dp * 10 ** self.log_Kpw)
                                             * self.permeation_surface_area 
                                             * delta_c / self.diffusion_path_length 
                                             * pipe.stagnation_time * self.stagnation_factor) 
                                            / self.ASSESSMENT_FACTOR_GROUNDWATER)

