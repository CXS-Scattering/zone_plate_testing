
# coding: utf-8

# In[1]:

import numpy as np
from skimage import io
import prop
import matplotlib.pyplot as plt
import pickle
import prop


# Importing all the required libraries. Numba is used to optimize functions.

# In[2]:

zp = np.load('zp.npy')
parameters = pickle.load(open('parameters.pickle','rb'))
parameters


# Importing zone plate pattern and the parameters associated with it.

# In[3]:

'''
decide whether to use TF or IR approach depending on the distance
'''
def decide(step_z,step_xy,L,wavel):
    dist = step_z
    sampling = step_xy
    critical = wavel*dist/L
    if sampling > critical :
        p = prop.propTF
        print('TF')
    else :
        p = prop.propIR
        print('IR')
    return p   
'''
used as part of the multislice loop
'''
def modify(wavefront,zp_delta,zp_beta,step_z,wavel):
    dist = step_z
    kz = 2 * np.pi * dist /wavel
    beta_slice = zp_beta
    delta_slice = zp_delta
    wavefront *= np.exp((kz * delta_slice) * 1j) * np.exp(-kz * beta_slice)
    return wavefront
'''
perform free space propogation using the method decided above
'''
def propogate(wavefront,step_xy,step_z,L,wavel,p):
    sampling = step_xy
    dist = step_z
    return p(wavefront,sampling,L,wavel,dist)


# *decide* : decides whether TF or IR approach should be used for propogation
# * *Inputs* : step size in z, step size in xy, support length, wavelength
# * *Outputs* : propogator
#     
# *modify* : wavefront is modified according to the material present
# * *Inputs* : wavefront, slice properties (here the zone plate),step size in z , wavelength
# * *Outputs* : modified wavefront
# 
# *propogate* : wavefront is propogated for the specified distance
# * *Inputs* : wavefront, step size in z, step size in xy, wavelength, propogator
# * *Outputs* : wavefront at output plane
# 

# In[4]:

zp_beta  = parameters['beta']*zp
zp_delta = parameters['delta']*zp
step_xy = parameters['step_xy']
energy = parameters['energy(in eV)']
wavel = parameters['wavelength in m']
wavefront = np.ones(np.shape(zp),dtype='complex64') #initialize wavefront 
prop_steps = 7500           #Total number of propogation steps
zone_plate_length = 1       #Number of steps along z for which zone plate is present
step_z = 1e-6               #Step size along z
L = step_xy*np.shape(zp)[0] #Support length
n = np.shape(zp)[0]         #dimension of input zp


# Setting up the parameters for the simulation

# In[7]:

'''
multislice simulation.
'''
def simulate(wavefront,zp,prop_steps,zone_plate_length,n,N1=0,N2=0):
        p = decide(step_z,step_xy,L,wavel)
        output = np.zeros((n,n,N2-N1))
        j = 0
        for i in range(prop_steps):
            print(i)
            if i<=zone_plate_length :
                wavefront = modify(wavefront,zp_delta,zp_beta,step_z,wavel)
                wavefront = propogate(wavefront,step_xy,step_z,L,wavel,p)
            else :
                wavefront = propogate(wavefront,step_xy,step_z,L,wavel,p)
            if i>N1 & i<N2:
                output[j] = abs(wavefront)
                j = j + 1
        return wavefront,output


# *simulate* : perform multislice simulation for the given zone plate
# * *Inputs* : wavefront, step size in z, step size in xy, wavelength, propogator
# * *Outputs* : wavefront at output plane
# 

# In[8]:

#simulate(wavefront,zp,prop_steps,zone_plate_length,n,5500,6500)

