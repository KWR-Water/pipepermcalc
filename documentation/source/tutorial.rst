========
Tutorial
========

Steps
-----

The following steps are needed to calculate the permeation of an organic chemical to a pipe using pipepermcalc:

#. Define segment(s) of a pipe using the segment class (pipe material, size).
#. Create a pipe from the segment(s) 
#. Set the groundwater contamination conditions (chemical, groundwater concentration and temperature) and mean daily flow rate.

Once a pipe has been created and the groundwater conditions and flow rate set the following calculations are possible:

#. Calculate the maximum or mean daily concentration in drinking water from a known groundwater concentration
#. Calculate the maximum or mean daily allowable concentration in groundwater to meet a drinking water norm for a given chemical

Example 1 - Simple
--------------------------------
Single pipe segment and chemical.

We always start by importing pipepermcalc Pipe() and Segment classes:

.. ipython:: python

    from pipepermcalc.pipe import * 
    from pipepermcalc.segment import * 

Step 1: Define the pipe segment(s) 
==================================
For this example there is only one pipe segment made of PE40. We define the length, inner diameter and thickness of the pipe in meters.

.. ipython:: python
    
    seg1 = Segment(name='seg1',
                    material='PE40',
                    length=25,
                    inner_diameter=0.0196,
                    thickness=0.0027,
                    )
Step 2: Create a pipe from the segment(s)
=========================================
We create a pipe from the segment using the Pipe() class by inputing the list of segment name(s).

.. ipython:: python

    pipe1 = Pipe(segment_list=[seg1])

Step 3: Set the groundwater conditions and flow rate
====================================================
Next we define the groundwater conditions. These are necessary if we want to calculate the concentration in drinking water for a given groundwater concentration.
We define the chemical of interest (name in Dutch), the concentration in groundwater (g/m3) and temperature (degrees Celcius). 
The name of the chemical is checked against the chemical database and the closest matching chemical is printed.
Finally we define the mean daily flow rate in m3/day.

.. ipython:: python
    
    pipe1.set_groundwater_conditions(chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8)

    pipe1.set_flow_rate(flow_rate=0.5)

Step 4: Calculate the drinking water concentration
==================================================
For the given groundwater conditions we can calculate the peak and mean daily concentration in drinking water for the pipe. 
The peak concentration is calculated as the concentration after a stagnation period (no flow in the pipe), with a default stagnation time of 8 hours is used.

.. ipython:: python
    
    peak_conc = pipe1.calculate_peak_dw_concentration(stagnation_time_hours = 8)
    print("The peak concentration is:", round(peak_conc,4), "g/m3")

    mean_conc = pipe1.calculate_mean_dw_concentration()
    print("The mean daily concentration is:", round(mean_conc,4), "g/m3")
                         
Step 5: Calculate the allowable groundwater concentration
=========================================================
It is also possible to calculate the allowable groundwater concentration which would not result in a concentration in drinking water exceeding a specified drinking water concentration for the given chemical.
Both the groundwater concentration which would not exceed the peak and the mean daily concentration can be calculated.
We define the chemical of interest, the target drinking water concentration (in this case the drinking water norm value), the temperature of the groundwater and, in the case of the peak concentration, the stagnation time.

.. ipython:: python

    peak_conc = pipe1.calculate_peak_allowable_gw_concentration(concentration_drinking_water=0.001,
                                stagnation_time_hours = 8,
                                chemical_name="Benzeen", 
                                temperature_groundwater=12)    
   
    print("The peak groundwater concentration, not exceeding the norm:", round(peak_conc,4), "g/m3")

    mean_conc = pipe1.calculate_mean_allowable_gw_concentration(concentration_drinking_water=0.001,
                                chemical_name="Benzeen", 
                                temperature_groundwater=12)    
   
    print("The mean groundwater concentration, not exceeding the norm:", round(mean_conc,4), "g/m3")


Miscellaneous
=============
The chemical/permeability information for the pipe can be inspected using the pipe_permeability_dict:

.. ipython:: python

    pipe1.pipe_permeability_dict

The individual segment information, e.g. volume, permeation surface area, logK, LogD etc., are attributes of the segments themselves:

.. ipython:: python

    seg1.volume

    seg1.permeation_surface_area

    seg1.log_Dp

    seg1.log_Kpw


Example 2 - Multiple segments
--------------------------------
In this example we create a pipe made from multiple segments with different permeation directions.

Depending on the types of pipe segment, the permeation direction can either be perpendicular (default) or parallel to the flow direction in the pipe. The diffusion path length is the length of permeation through the pipe segment.

.. image:: images/pipe_schematic.png
  :width: 600
  :alt: pipe_schematic.png

In scenarios 1 and 3 above, the permeaiton is perpendicular to the flow direction and the volume is calculated from the segment dimensions. The surface area is given as the inner surface area of the segment. In pipepermcalc the default permeation direction is perpendicular and the diffusion path length equal to the thickness of the pipe length.

In the example shown above, permeation is *parallel* to the flow direction through a connecting rubber in scenario 2. For this scenario, the volume is assumed to be zero and the permeation surface area is the annular area of the rubber. The diffusion path length in this case is equal to the length of the segment.

In the following example we create a pipe made from two 5m PE40 pipe segments, joined by a EPDM ring with permeation parallel to the flow direciton:

.. ipython:: python

    seg1 = Segment(name='seg1',
                material='PE40',
                length=5,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

    seg2 = Segment(name='seg2',
                    material = 'EPDM',
                    length=0.06,
                    inner_diameter=0.025,
                    thickness=0.001,
                    diffusion_path_length = 0.06, 
                    permeation_direction = 'parallel'
                    )

    seg3 = Segment(name='seg3',
                material='PE40',
                length=5,
                inner_diameter=0.0196,
                thickness=0.0027,
                )

    pipe2 = Pipe(segment_list=[seg1, seg2, seg3])


As seen in the example above, only the segment with the parallel flow requires a specified permeation direction, as the default is perpendicular, and the diffusion path length, as the default is the thickness.

The remaining calculations are done the same as for the simple example:

.. ipython:: python

    pipe2.set_groundwater_conditions(chemical_name="Benzeen", 
                                temperature_groundwater=12, 
                                concentration_groundwater=1.8)

    pipe2.set_flow_rate(flow_rate=0.5)

    peak_conc = pipe2.calculate_peak_dw_concentration(stagnation_time_hours = 8)
    print("The peak concentration is:", round(peak_conc,4), "g/m3")

    mean_conc = pipe2.calculate_mean_dw_concentration()
    print("The mean daily concentration is:", round(mean_conc,4), "g/m3")


Example 3 - Specify the K and D used
------------------------------------



.. Advanced tutorials: 
.. 1) set tolerance, relaxation factor, max iterations etc., 
.. 2) mutliple segments w/diagram (also add to glossary?), 
.. 3) change the K, D used
.. other: view the pipe segment information, view pipe permeability information 