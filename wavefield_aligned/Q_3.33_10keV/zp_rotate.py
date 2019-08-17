# Script to perform rotation of zone plate array.
# Expand the number of slices along z direction so that
# step size along axis of rotation matches step size
# along direction of propagation. Then use pyvips to 
# perform the rotation and finally collapse down to
# the number of slices as required.

import numpy as np
import matplotlib.pyplot as plt
import skimage.transform
import h5py,sys
from tqdm import tqdm_notebook,trange
from joblib import Parallel, delayed

# function to fill the expanded array with
# un-tilted zone plate.
def expand(j,loc,A,B):
    A[:,:,loc-j] = B
    A[:,:,loc+j] = B
    return


# This file accepts one input parameter : the angle
# by which the zone plate must be tilted by. Everything
# else is read from the hdf5 dataset attributes! 

if __name__=="__main__":
    
    # Read the input arguments and set rotation angle. 
    angle = float(sys.argv[1])*0.1
    print(angle)

    # Load the pyvips conversion dictionary.
    get_ipython().run_line_magic('run', 'pyvips_dict.py')

    # Load the reduced zone plate pattern.
    f = h5py.File("zp_reduced.hdf5","r")
    zp_dset = f['zp_reduced']
    zp_reduced = zp_dset[:,:]
    
    # Load the parameters associated with the zone 
    # plate from the attributes.
    suggested_step_size = zp_dset.attrs['suggested_step_size']
    zp_thickness        = zp_dset.attrs['zp_thickness']
    rot_slices          = zp_dset.attrs['rotation_slices']
    step_xy             = zp_dset.attrs['step_xy']
    expand_slices       = zp_dset.attrs['expand_slices']
    step_z              = zp_dset.attrs['step_z']
    max_angle           = zp_dset.attrs['max_angle']
    
    # Close the file.
    f.close()

    # Print the parameters.
    print("rotation_slices     :",rot_slices)
    print("step_xy             :",step_xy)
    print("expand_slices       :",expand_slices)
    print("step_z              :",step_z)
    print("zp_thickness        :",zp_thickness)

    
    # Remove old copies and create a new hdf5 file and dataset
    # to hold the final (sub) array of refractive indices. Ensure
    # that the dataset is chunked for fast I/O along both the 
    # axis of rotation and axis of beam propagation.
    get_ipython().system('rm zp_rotate.hdf5')
    f = h5py.File("zp_rotate.hdf5")
    
    zp_rotate = f.create_dataset("zp_rotate",
                            (np.shape(zp_reduced)[0],np.shape(zp_reduced)[1],rot_slices),
                            chunks=True)

    # Store the parameters as dset attributes.
    zp_rotate.attrs['suggested_step_size'] = suggested_step_size
    zp_rotate.attrs['zp_thickness']        = zp_thickness
    zp_rotate.attrs['step_xy']             = step_xy
    zp_rotate.attrs['rotation_slices']     = rot_slices
    zp_rotate.attrs['expand_slices']       = expand_slices
    zp_rotate.attrs['step_z']              = step_z
    zp_rotate.attrs['angle']               = angle


    # Calculate the number of slices the zone plate
    # would occupy in the expanded array.
    # Convention : zone plate would occupy 2x zp_extent
    # slices, centered around center_loc.
    zp_extent  = int((zp_thickness/step_xy)/2)
    center_loc = int(expand_slices*rot_slices/2)

    # Start the main rotation loop. We will rotate 500 slices
    # along the axis of rotation at one time. Each rotate iteration
    # is expand -> rotate -> collapse. Expand and collapse use
    # standard NumPy functions. Rotate uses the python binding of 
    # libvips (and can rotate a stack of 2D images with multiple 
    # threads).

    for i in tqdm_notebook(range(30)):

        # Create temporary arrays to hold the compact and 
        # expanded pattern. Create a transfer array to 
        # hold the set of zone plate indices currently being 
        # rotated.
        temp  = np.zeros((500,15000,rot_slices))
        temp_ = np.zeros((500,15000,expand_slices*rot_slices))
        transfer = np.zeros((500,15000))

        # Hold no tilt pattern for initial transfer
        transfer = zp_reduced[500*i:500*(i+1),:]

        #Fill expanded array
        Parallel(n_jobs=24,backend='threading',require='sharedmem')\
        (delayed(expand)(j,center_loc,temp_,transfer) for j in trange(zp_extent,desc='Expand   '))
        
        # Rotate expanded array when angle is set to 
        # anything other than zero. Use pyvips to perform the
        # rotation, which automatically uses all available
        # threads.
        if angle!=0:
            for j in trange(500,desc='Rotate zp'):
                temp_[j,:,:] = pyvips_rotate(temp_[j,:,:],angle_ = angle)

        # Collapse down to number of number of slices
        # that will be used for beam propagation.
        for j in trange(rot_slices,desc='Collapse '):
            temp[:,:,j] = np.sum(temp_[:,:,expand_slices*j:expand_slices*(j+1)],axis=-1)/expand_slices

        # Transfer the collapsed pattern to the hdf5 dset
        zp_rotate.write_direct(temp, np.s_[:,:,:], np.s_[500*i:500*(i+1),:,:])

        # Delete all temporary arrays.
        del temp_
        del temp
        del transfer

    # After rotation is complete, set the angle 
    # attribute and close the file.
    zp_rotate.attrs['angle'] = angle
    f.close()
