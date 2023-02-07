========
Tutorial
========
Practical examples


We always start by importing package_name:

.. ipython:: python

    import package_name

Steps
-----

The following steps are needed to calculate the permeation of an organic chemical to a pipe using Pipepermcalc:

#. Define a pipe using the Pipe class. 
#. Set the groundwater contamination conditions (chemical, groundwater concentration and temperature)
#. Define the different segments of the pipe of interest (pipe material, size)
#. Calculate the partitioning and diffusion coefficients for the pipe material and chemical of interest
#. Calculate the maximum or mean daily concentration in the groundwater and soil

Example 1
--------------------------------

Single pipe segment and chemical.

.. ipython:: python

   from pipepermcalc.ppc_overview_script import * 

Step 1: Define a Pipe and set groundwater conditions
==============================================

.. ipython:: python
    
    pipe1 = Pipe()
    pipe1.set_groundwater_conditions(chemical_name="Benzene", 
                                 temperature_groundwater=12, 
                                 concentration_groundwater = 1.8)

The units for the concentration in groundwater are g/m3 and temperature are degrees Celcius.
The chemicals available are XX
 .. @Bram something here about the enum chemical name search

Step 2: Define the pipe segments
==============================================
For this example there is only one pipe segment made of PE40. We define the length, diameter and thickness of the pipe in meters and the flow rate through the pipe in m3/day.

.. ipython:: python

    pipe1.add_segment(name='segment_1',
                    material='PE40',
                    length=25,
                    diameter=0.0196,
                    thickness=0.0027,
                    flow_rate=0.5)

To view the pipe information and different pipe segments:

.. ipython:: python
    pipe1.pipe_dictionary


Step 3: Calculate the partitioning and diffusion coefficients
==============================================
(optional?) step to calculate the partitioning and diffusion coefficients for the pipe segment. These depend on both the chemical of interest, pipe material.
Corrections for the groundwater concentration and temperature are taken into account.

.. ipython:: python

    pipe1.calculate_pipe_K_D(
                    pipe_material= "PE40")

.. @ah_todo update this is we move to calculating per segment instead of defining the pipe material

To view the calculated partitioning and diffusion coefficients:

.. ipython:: python
    pipe1.pipe_permeability_dict


Step 4: Calculate the peak drinking water concentration
==============================================
To calculate the peak drinking water concentration in a pipe segment, a stagnation period is defined. This is a period of time, in hours, in which the water in the pipe is not flowing, for example when no water is used at night. A default stagnation period of 8 hours is used.

.. ipython:: python

    pipe1.calculate_max_dw_concentration(stagnation_time_hours = 8, 
                                    pipe_segment='seg1')


To view the peak drinking water concentration after the stagnation period:

.. ipython:: python

    pipe1.pipe_permeability_dict['concentration_peak_after_stagnation']