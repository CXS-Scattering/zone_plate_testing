#!/bin/bash
#SBATCH --job-name=zone_plate_testing
#SBATCH --nodes=1 
#SBATCH -p apsxrmd
#SBATCH --time=120:00:00

~/.conda/envs/intelpy3/bin/python simulate_var_width_thickness.py
