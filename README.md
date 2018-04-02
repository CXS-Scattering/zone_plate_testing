# zone_plate_testing

### Workflow for creating Fresnel zone plate patterns and studying their performance. The input wavefront can be tilted to simulate the effect of zone plate misalignment.

*create_rings.ipynb* : creates rings according to the input parameters. 

*simulate_zp_with_tilt.ipynb* : performs multislice simulation to evaluate the performance of the zone plate with the option to tilt the input wave. 

*simulate_var_width_thickness.py* : python script to test effects of tile on zone plate performance for a various thicknesses and number of zones. Essentially an extenstion of the jupyter notebook above. 

*analyze_tilt.ipynb*              : analyze the results from running the above script. 

### Example
The results from a test run for a zone plate with outermost zone width : 24 nm, radius : 12 um for a thickness of 2 um at 10 keV are shown.

### Dependencies 
Uses the [multislice](https://github.com/sajid-ali-nu/multislice/) module.
