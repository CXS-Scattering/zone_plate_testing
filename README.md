## Zone Plate Testing

#### Currently refactoring code. Check back in a few weeks.

Workflow for creating Fresnel zone plate patterns and studying their performance. Simulations of zone plate tilt misalignment.

`create_rings.ipynb` : creates rings according to the input parameters.

`create_rings.py` : python script that runs the same as above.

`simulate_zp_with_tilt.ipynb` : performs multislice simulation to evaluate the performance of the zone plate with the option to tilt the input wave. 

`simulate_var_width_thickness.py` : python script to test effects of tile on zone plate performance for a various thicknesses and number of zones. Essentially an extenstion of the jupyter notebook above. 

`analyze_tilt.ipynb`              : analyze the results from running the above script. 

### Caution
The code retreives the refractive indices only for elements and not compounds. If one intends to use compounds then the delta and beta values need to be added manually. 

### Acknowledgements
“We gratefully acknowledge the computing resources provided on Blues (and/or Bebop), a high-performance computing cluster operated by the Laboratory Computing Resource Center at Argonne National Laboratory.”
