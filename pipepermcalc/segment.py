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
        
