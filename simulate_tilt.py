import numpy as np
import pickle
from skimage import io,img_as_float
import time
from multislice import prop,prop_utils

def find_edge(x):
    if x<7500:
        if x>100:
            return 100
        else :
            return int(np.floor(x / 2) * 2)
    else :
        if (15000-x)>100:
            return 100
        else :
            return int(np.floor(15000-x / 2) * 2)

def get_focal_spot(focal_plane_):
    x_,y_ = np.where(focal_plane_==np.max(focal_plane_))
    x_ = x_[0]
    y_ = y_[0]
    print(x_,y_)
    x1 = find_edge(x_)
    y1 = find_edge(y_)
    focal_spot_ = np.zeros((200,200))
    if (x1+y1) != 200:
        focal_spot_[100-x1:100+x1,100-y1:100+y1] = focal_plane_[x_-x1:x_+x1,y_-y1:y_+y1]
    else :
        focal_spot_[:,:] = focal_plane_[x_-100:x_+100,y_-100:y_+100]
    return focal_spot_,x_,y_,np.max(focal_plane_)


def tilt(i):
    t1 = time.time()
    print('claclulating for tilt angle : ',i)      
    theta = (i)*(np.pi/180)
    slope = np.tan(theta)
    x = np.linspace(zp_coords[0]*1e-6,zp_coords[1]*1e-6,n)
    X,Y = np.meshgrid(x,x)
    z1 = 2*np.pi*(1/wavel)*slope*X
    wave_in = np.multiply(np.ones((n,n),dtype='complex64'),np.exp(1j*(z1)))
    number_of_steps_zp =  prop_utils.number_of_steps(step_xy,wavel,zp_thickness) + 1
    wave_focus = prop_utils.propogate_through_object(wave_in,zp_delta,zp_beta,zp_thickness,step_xy,wavel,L,
                                                     number_of_steps_zp,d1=0,d2=f, xray_object='zone Plate')
    focal_spot,x_,y_,max_val = get_focal_spot(np.abs(wave_focus))
    io.imsave('tilt_'+str(i)+'.tiff',img_as_float(focal_spot))
    t2 = time.time()
    print('tilt image number :',i,'time taken :',(t2-t1))
    print(x_,y_,max_val)
    return {'tilt_angle':i,'max_loc':np.array([x_,y_]),'max_val':max_val}
    
zp = np.load('Au_417_zp.npy')
parameters = pickle.load(open('parameters_Au_417.pickle','rb'))
zp_coords = parameters['zp_coords']    
zp_beta  = parameters['beta']*zp
zp_delta = parameters['delta']*zp
step_xy = parameters['step_xy']
energy = parameters['energy(in eV)']
wavel = parameters['wavelength in m']
f = parameters['focal_length']
L = step_xy*np.shape(zp)[0] 
n = np.shape(zp)[0]         
zp_thickness = 1e-6    
print('all imports done')
time.sleep(1)
inputs = np.linspace(-1,1,50)
max_loc = []
for i in inputs :
    max_loc.append(tilt(angle,zp,t,parameters,out_dir))
np.save('max_loc.npy', max_loc) 
