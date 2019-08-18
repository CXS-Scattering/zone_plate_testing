# Script to perform multislice simulation
# after the rotation has been performed.
# Read the (sub) pattern from hdf5 file and 
# fill in the empty space with zeros. 
# Perform multi-slice propagation to get exit
# wave. Propagate to focus - offset. Offset is 
# half the thickness of the zone plate. 
# Save result as .npy files.

import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
import pickle,os
from multislice import prop,prop_utils
from os.path import dirname as up
from tqdm import tqdm_notebook
import h5py,sys

# This file accepts one input parameter : the angle
# by which the zone plate must be tilted by. Everything
# else is read from the hdf5 dataset attributes! 

if __name__=="__main__":

    # Read the input arguments and set rotation angle. 
    angle = float(sys.argv[1])*0.1
    print(angle)

    # Read the relevant parameters saved when the 
    # zone plate pattern was created. Switch back 
    # to current directory.
    pwd = os.getcwd()
    os.chdir(up(up(os.getcwd()))+str('/zp_database/hard_xray_zp'))
    
    parameters = pickle.load(open('parameters.pickle','rb'))
    for i in parameters : print(i,' : ',parameters[i])
    
    wavel = parameters['wavelength in m']
    step_xy = parameters['step_xy']
    zp_coords = parameters['zp_coords']
    grid_size = parameters['grid_size']
    beta  = parameters['beta']
    delta = parameters['delta']
    step_xy = parameters['step_xy']
    wavel = parameters['wavelength in m']
    f = parameters['focal_length']
    L = step_xy*grid_size
    
    os.chdir(pwd)
    
    # This functions reads the hdf5 to extract
    # the rotated zone plate pattern at current
    # slice index. Fill empty space with zeros.
    # Will be called at each slice of the
    # multi-slice loop.

    def get_pattern(i):
        temp  = np.zeros((55296,55296))
        temp_ = np.zeros((15000,15000))
        f = h5py.File("zp_rotate.hdf5")
        zp_rotate = f['zp_rotate']
        zp_rotate.read_direct(temp_, np.s_[:,:,i], np.s_[:,:])
        f.close()
        N = int(55296/2)
        n = 7500
        temp[N-n:N+n,N-n:N+n] = temp_
        del temp_
        return temp
    
    # Create the input wave.
    wavefront = np.ones((55296,55296),dtype='complex64')

    # Read additional parameters from hdf5 file.
    f = h5py.File("zp_rotate.hdf5")
    zp_rotate = f['zp_rotate']
    number_of_steps = zp_rotate.attrs['rotation_slices']
    step_z = zp_rotate.attrs['step_z']
    f.close()

    # Set number of threads for numexpr.
    ne.set_vml_num_threads(24)

    # Pick the appropriate propagator. 
    p = prop_utils.decide(step_z,step_xy,L,wavel)

    # Main multi-slice loop. Extract pattern at 
    # each slice. Propagate until exit plane.
    for i in tqdm_notebook(range(number_of_steps)):
        pattern = get_pattern(i)
        wavefront = prop_utils.modify_two_materials_case_2(wavefront,step_z,wavel,
                                                pattern,delta,beta,
                                                np.ones(np.shape(pattern))-pattern,0,0)
        wavefront,L1 = p(wavefront,step_xy,L,wavel,step_z)

    wave_exit = wavefront
    del wavefront

    # Propagate to focus.
    step_z = f - (number_of_steps*step_z)/2
    p = prop_utils.decide(step_z,step_xy,L,wavel)
    print('Propagation to focal plane')
    print('Fresnel Number :',((L**2)/(wavel*step_z)))
    wave_focus,L2 = p(wave_exit - np.ones(np.shape(wave_exit)),step_xy,L,wavel,step_z)
    wave_focus = wave_focus +  np.ones(np.shape(wave_exit))

    # Extract the focal spot and save it.
    focal_spot_size = 250
    focal_spot,x_,y_,max_val = prop_utils.get_focal_spot(wave_focus,grid_size,focal_spot_size)

    np.save('foc_spot_Q_3.33_'+str(round(angle,3))+'_degree.npy',focal_spot)
    np.save('foc_loc_Q_3.33_'+str(round(angle,3))+'_degree.npy',np.array([x_,y_]))
    



