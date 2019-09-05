'''

This file is part of the zone_plate_testing repository and an exntension of the recipie to simulate one tilted zone plate fonund in the notebook simlate_zp_with_tilt.ipynb to a workflow which does the aforementioned task multiple times. 

'''

import numpy as np
import numexpr as ne
import os,pickle,time
import matplotlib.pyplot as plt
from skimage import io
from skimage import img_as_float
from multislice import prop,prop_utils,fft_utils
from os.path import dirname as up
from tqdm import tqdm, trange

def propTF_modified(u,step,L1,wavel,z,n,FX,FY,fft_object = None) :
    pi = np.pi
    if fft_object != None :
        fft_object.run_fft2(u)
    else: 
        u = np.fft.fft2(u)
    
    _wavel = wavel*(1/n)
    
    u = ne.evaluate('exp(-1j*(2*pi*z/_wavel)*sqrt(1-_wavel**2*(FX**2+FY**2)))*u')
    
    if fft_object != None :
        fft_object.run_ifft2(u)
    else :
        u = np.fft.ifft2(u)
    
    return u,L1


# In[3]:


def modify_two_materials_case_2_new(wavefront,step_z,wavel,pattern_1,delta_1,beta_1,pattern_2,delta_2,beta_2):
    pi = np.pi
    kz = ne.evaluate('2 * pi * step_z /wavel')
    
    #Collapsing the following into one numexpr statement : 
    modulation_1 = ne.evaluate('exp(kz*((1 - delta_1 + 1j*beta_1)**2 - 1))')
    modulation_2 = ne.evaluate('exp(kz*((1 - delta_2 + 1j*beta_2)**2 - 1))')
    return ne.evaluate('wavefront * ( pattern_1*modulation_1 + pattern_2*modulation_2 )')
    
    
    #return ne.evaluate('wavefront * ( pattern_1*exp((kz*delta_1)*1j - kz*beta_1)+pattern_2*exp((kz*delta_2)*1j - kz*beta_2) )')


# In[4]:


def optic_illumination_modified(wavefront_input,
                       pattern,delta,beta,
                       thickness,step_xy,wavel,
                       number_of_steps,d1,d2,n,use_fftw='False',**kwargs):
    
    wavefront = np.copy(wavefront_input)
    L = np.shape(wavefront_input)[0]*step_xy
    xray_object = str('zone plate')
    mode = str('serial')
    if use_fftw == 'True':
        fft_obj = fft_utils.FFT_2d_Obj(np.shape(wavefront),flag='PATIENT')
    else :
        fft_obj = None

    if 'xray_object' in kwargs :
        xray_object = kwargs['xray_object']
    if 'mode' in kwargs : 
        mode = kwargs['mode']    
      
    
    #pre object
    if d1 != 0 :
        print('Free space propogation before '+str(xray_object)+'...')
        step_z = d1
        p = prop_utils.decide(step_z,step_xy,L,wavel)
        print('Fresnel Number :',((L**2)/(wavel*step_z)))
        wavefront,L  = p(wavefront,step_xy,L,wavel,step_z,fft_obj)
    
    #through object
    step_z = thickness/number_of_steps
    for i in tqdm(range(number_of_steps),desc='Propogation through '+str(xray_object)+'...'):
        wavefront = modify_two_materials_case_2_new(wavefront,step_z,wavel,pattern,delta,beta,np.ones(np.shape(pattern))-pattern,0,0)
        M,N = np.shape(pattern)
        fx,fy = np.meshgrid(np.fft.fftfreq(M,step_xy),np.fft.fftfreq(N,step_xy))
        wavefront,L  = propTF_modified(wavefront,step_xy,L,wavel,step_z,1,fx,fy)
        del fx,fy

    #post object
    if d2 !=0 :
        step_z = d2
        print('Free space propogation after '+str(xray_object)+'...')
        p = prop_utils.decide(step_z,step_xy,L,wavel)
        print('Fresnel Number :',((L**2)/(wavel*step_z)))
        wavefront,L  = p(wavefront,step_xy,L,wavel,step_z)
    
    wavefront_out = np.copy(wavefront)
    del wavefront
    return wavefront_out,L

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

tilt : get the focal spot for a given zone plate and tilt angle of the input wave and save it to a tiff file

Inputs : i-tilt angle in degrees,zp - zone plate pattern, thickness (of the zone plate),parameters, out_dir - output directory

Outputs : dictionary containing :  tilt angle, max_loc - location of the zone plate, L2 - side length at output plane

'''

def tilt(i,zp,thickness,parameters,out_dir):
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
    
    os.chdir(out_dir)

    t1 = time.time()
    print('claclulating for tilt angle : ',i)      
    theta = (i)*(np.pi/180)
    slope = np.tan(theta)
    x = np.linspace(zp_coords[0],zp_coords[1],n)
    X,Y = np.meshgrid(x,x)
    z1 = 2*np.pi*(1/wavel)*slope*X
    wave_in = np.multiply(np.ones((n,n),dtype='complex64'),np.exp(1j*(z1)))
    ne.set_vml_num_threads(24)
    print('step_xy :',step_xy,'wavelength :',wavel,'zp_thickness :',zp_thickness)
    number_of_steps_zp =  (prop_utils.number_of_steps(step_xy,wavel,zp_thickness))*2
    wave_exit,L2 = optic_illumination_modified(wave_in,zp,delta,beta,zp_thickness,step_xy,wavel,number_of_steps_zp,0,0,1-delta+1j*beta)
    step_z = f
    p = prop_utils.decide(step_z,step_xy,L,wavel)
    print('Propagation to focal plane')
    print('Fresnel Number :',((L**2)/(wavel*step_z)))
    wave_focus,L2 = p(wave_exit - np.ones(np.shape(wave_exit)),step_xy,L,wavel,step_z)
    #wave_focus,L2 = prop.prop1FT(wave_exit - np.ones(np.shape(wave_exit)),step_xy,L,wavel,step_z)
    #wave_focus = wave_focus +  np.ones(np.shape(wave_exit))
    focal_spot,x_,y_,max_val = prop_utils.get_focal_spot(np.abs(wave_focus),grid_size)
    focal_spot = np.abs(focal_spot)
    io.imsave('tilt_'+str(i)+'.tiff',img_as_float(focal_spot))
    t2 = time.time()
    
    print('tilt image number :',i,'time taken :',(t2-t1))
    
    return {'tilt_angle':i,'max_loc':np.array([x_,y_]),'L2':L2}


os.chdir(os.getcwd()+str('/rings'))
parameters = pickle.load(open('parameters.pickle','rb'))
grid_size = parameters['grid_size']

'''
 Creating the set of parameters for the simulation data set. Making a common output directory for the results.
'''

num_zones = np.array([500])
##thickness = np.array([3.4227e-08]) #Q=0.333..
##thickness = np.array([0.342269e-6])  #Q=3.333..
thickness = np.array([1.0268e-6]) #Q=10..
inputs = np.linspace(0,5,51)


output_dir = up(os.getcwd())+str('/output/')
try:
    os.chdir(output_dir)
except:
    os.mkdir(output_dir)

'''
Nested for loop to run the simluation over the variable number of zones, thicknesses and tilt angles. The output for each set of parameters is saved into a different folder.
'''

for var1 in range(len(num_zones)):
    n = num_zones[var1]
    os.chdir(str('/home/sajid/Documents/ZP_simulations/var_width_thickness_final/0.5kev/rings'))
    zp =  make_zp_from_rings(n,grid_size)
    os.chdir(str('/home/sajid/Documents/ZP_simulations/var_width_thickness_final/0.5kev/output'))
    for var2 in range(len(thickness)) :
        print(num_zones[var1],' zones, ',thickness[var2] * 1e6,' microns thick')
        max_loc = []
        t = thickness[var2]
        out_dir = os.getcwd()+ str('/zones_'+str(n)+'_thickness_'+str(t))
        os.mkdir(out_dir)
        print(n,t,grid_size)
        for angle in inputs:
            max_loc.append(tilt(angle,zp,t,parameters,out_dir))
        np.save('max_loc.npy', max_loc) 
