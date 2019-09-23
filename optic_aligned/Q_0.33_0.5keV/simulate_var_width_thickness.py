'''

This file is part of the zone_plate_testing repository and an exntension of the recipie to simulate one tilted zone plate fonund in the notebook simlate_zp_with_tilt.ipynb to a workflow which does the aforementioned task multiple times. 

'''

import numpy as np
import os,pickle
import matplotlib.pyplot as plt
from multislice import prop,prop_utils
from os.path import dirname as up


'''

make_zp_from_rings : make a zone plate from the rings which were created earlier.

Inputs : n - number of rings, grid_size 

Outputs : a numpy array containing the zone plate

'''
def make_zp_from_rings(n,grid_size):
    zp = np.zeros((grid_size,grid_size))
    for i in range(n):
        if i%2 == 1 :
            locs_ = np.load('ring_locs_'+str(i)+'.npy')
            locs_ = tuple((locs_[0],locs_[1]))
            vals_ = np.load('ring_vals_'+str(i)+'.npy')
            zp[locs_] = vals_
    return zp


'''

tilt : get the focal spot for a given zone plate and tilt angle of the input wave and save it.

Inputs : i-tilt angle in degrees,zp - zone plate pattern, thickness (of the zone plate), parameters.

'''

def tilt(i,zp,thickness,parameters):
    zp_thickness = thickness 
    beta  = parameters['beta']
    delta = parameters['delta']    
    zp_coords = parameters['zp_coords']    
    step_xy = parameters['step_xy']
    energy = parameters['energy(in eV)']
    wavel = parameters['wavelength in m']
    f = parameters['focal_length']
    L = step_xy*np.shape(zp)[0] 
    n = np.shape(zp)[0]   

    print('claclulating for tilt angle : ',i)      
    theta = (i)*(np.pi/180)
    slope = np.tan(theta)
    x = np.linspace(zp_coords[0],zp_coords[1],n)
    X,Y = np.meshgrid(x,x)
    z1 = 2*np.pi*(1/wavel)*slope*X
    wave_in = np.multiply(np.ones((n,n),dtype='complex64'),np.exp(1j*(z1)))
    
    number_of_steps_zp =  (prop_utils.number_of_steps(step_xy,wavel,zp_thickness)+1)*2
    wave_focus,L2 = prop_utils.optic_illumination(wave_in,zp,delta,beta,zp_thickness,step_xy,wavel,number_of_steps_zp,0,f)
    focal_spot,x_,y_,max_val = prop_utils.get_focal_spot(wave_focus,grid_size)
    
    np.save('foc_spot_Q_0.33_'+str(round(angle,3))+'_degree.npy',focal_spot)
    np.save('foc_loc_Q_0.33_'+str(round(angle,3))+'_degree.npy',np.array([x_,y_]))
    
    return


'''
Load the zone plate from memeory, set parameters.
'''
pwd = os.getcwd()
os.chdir(up(up(os.getcwd()))+str('/zp_database/soft_xray_zp/'))
parameters = pickle.load(open('parameters.pickle','rb'))
grid_size = parameters['grid_size']
num_zones = 700
zp =  make_zp_from_rings(num_zones,grid_size)
os.chdir(pwd)

thickness = 0.03807e-6 #Q=0.333
inputs = np.arange(51)*0.1


print(num_zones,' zones, ',thickness*1e6,' microns thick')
max_loc = []

'''
Loop to run the simluation over the variable number of tilt angles.
'''
for angle in inputs:
    tilt(angle,zp,thickness,parameters)
