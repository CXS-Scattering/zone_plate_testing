#!/usr/bin/env python
import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
import pickle,os
from multislice import prop,prop_utils
from os.path import dirname as up
from tqdm import tqdm_notebook
import h5py,sys

if __name__=="__main__":

    angle = float(sys.argv[1])*0.1
    print(angle)

    pwd = os.getcwd()
    os.chdir(up(os.getcwd())+str('/rings'))
    parameters = pickle.load(open('parameters.pickle','rb'))
    for i in parameters : print(i,' : ',parameters[i])
    wavel = parameters['wavelength in m']
    step_xy = parameters['step_xy']
    zp_coords = parameters['zp_coords']
    grid_size = parameters['grid_size']
    os.chdir(pwd)


    beta  = parameters['beta']
    delta = parameters['delta']
    step_xy = parameters['step_xy']
    wavel = parameters['wavelength in m']
    f = parameters['focal_length']
    L = step_xy*55296 

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

    wavefront = np.ones((55296,55296),dtype='complex64')

    number_of_steps = 13
    step_z =  6.419657057352801e-07
    ne.set_vml_num_threads(12)

    p = prop_utils.decide(step_z,step_xy,L,wavel)

    for i in tqdm_notebook(range(number_of_steps)):
        pattern = get_pattern(i)
        #print("pattern extracted!")
        wavefront = prop_utils.modify_two_materials_case_2(wavefront,step_z,wavel,
                                                pattern,delta,beta,
                                                np.ones(np.shape(pattern))-pattern,0,0)
        wavefront,L1 = p(wavefront,step_xy,L,wavel,step_z)

    wave_exit = wavefront
    del wavefront

    step_z = f - (number_of_steps*step_z)/2
    p = prop_utils.decide(step_z,step_xy,L,wavel)
    print('Propagation to focal plane')
    print('Fresnel Number :',((L**2)/(wavel*step_z)))
    wave_focus,L2 = p(wave_exit - np.ones(np.shape(wave_exit)),step_xy,L,wavel,step_z)
    wave_focus = wave_focus +  np.ones(np.shape(wave_exit))

    focal_spot_size = 100
    focal_spot,x_,y_,max_val = prop_utils.get_focal_spot(wave_focus,grid_size,focal_spot_size)

    np.save('foc_spot_Q_3.33_'+str(round(angle,3))+'_degree.npy',focal_spot)
    



