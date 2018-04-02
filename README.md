## zone_plate_testing

Codes for creating Fresnel zone plate patterns and to study their performance. Effects of tilting the zone plate can also be studied.

*create_rings.ipynb*              - creates rings according to the input parameters. <br>
*simulate_zp_with_tilt.ipynb*     - performs multislice simulation to evaluate the performance of the zone plate with the option to tilt the input wave. <br>
*simulate_var_width_thickness.py* - python script to test effects of tile on zone plate performance for a various thicknesses and number of zones. Essentially an extenstion of the jupyter notebook above. <br>
*analyze_tilt.ipynb*              - analyze the results from running the above script. <br>

Dependencies :  uses the [multislice](https://github.com/sajid-ali-nu/multislice/) module.
