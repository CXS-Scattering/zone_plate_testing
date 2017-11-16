import numpy as np
import pickle
from skimage import io,img_as_float
import time
from multislice import prop,prop_utils

def get_focal_spot(focal_plane_):
    focal_spot = np.zeros((200,200))
    x_,y_ = np.where(focal_plane_==np.max(focal_plane_))
    x_ = x_[0]
    y_ = y_[0]
    if (x_ > 100) & (y_>100):
        n = 100
        focal_spot = focal_plane_[x_-n:x_+n,y_-n:y_+n]
    else : 
        n = 50
        focal_spot[50:150,50:150] = focal_plane_[x_-n:x_+n,y_-n:y_+n]
    return focal_spot


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
    focal_spot = get_focal_spot(np.abs(wave_focus))
    io.imsave('focal_spot_tile_anlge'+str(i)+'.tiff',img_as_float(focal_spot)
    t2 = time.time()
    print('tilt image number :',i,'time taken :',(t2-t1))


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
for i in inputs :
    tilt(i)

