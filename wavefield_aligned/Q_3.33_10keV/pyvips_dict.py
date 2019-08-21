# 2D/3D rotate utiliy function and 
# conversion dictionary.

import numpy as np
import pyvips

# conversion between numpy/pyvips.

format_to_dtype = {
    'uchar': np.uint8,
    'char': np.int8,
    'ushort': np.uint16,
    'short': np.int16,
    'uint': np.uint32,
    'int': np.int32,
    'float': np.float32,
    'double': np.float64,
    'complex': np.complex64,
    'dpcomplex': np.complex128,
}

dtype_to_format = {
    'uint8': 'uchar',
    'int8': 'char',
    'uint16': 'ushort',
    'int16': 'short',
    'uint32': 'uint',
    'int32': 'int',
    'float32': 'float',
    'float64': 'double',
    'complex64': 'complex',
    'complex128': 'dpcomplex',
}

# Utility function to convert 2D/3D
# numpy arrays using pyvips. 

def pyvips_rotate(x,angle_):
    # convert angle to radians
    # get cos/sin 
    angle_rad = angle_*(np.pi/180)
    tc = np.cos(angle_rad)
    ts = np.sin(angle_rad)
    
    # Set rows and columns for 2D transform
    # This would still work for a stack of 
    # 2D images
    rows, cols = x.shape[0], x.shape[1]

    # Set center of array to be the 
    # axis of rotation.
    c_ = np.array((rows,cols))/2 - 0.5

    # Check if we are working with a 
    # single 2D array or a stack of them
    if len(x.shape)!=3 : 
        height, width = x.shape
        bands = 1
    else :
        height, width,bands = x.shape
    
    # Create a pyvips image from the numpy array
    im = pyvips.Image.new_from_memory(x.reshape(width * height * bands).data, width, height, bands,
                                      dtype_to_format[str(x.dtype)])

    # Specify the interpolation type to bilinear. 
    # this interpolation will be used for roataion.
    inter = pyvips.vinterpolate.Interpolate.new('bilinear')

    # Perform the rotation via the pyvips
    # affine transform. Use the aforementioned 
    # bilinear interpolator.
    im = im.affine([tc, ts, -ts, tc], interpolate = inter,
                   idx=-c_[1],idy=-c_[0],
                   odx=c_[1],ody=c_[0],
                   oarea=[0, 0, im.width, im.height])

    # Transfer the pyvips image to numpy
    b = np.ndarray(buffer=im.write_to_memory(),dtype=format_to_dtype[im.format],shape=[im.height, im.width, im.bands])
    
    # delete the pyvips image 
    del im

    # Reshape the data to be either a 2D
    # 3D array.
    if len(x.shape)==3 : 
        b = b.reshape(np.shape(b)[0],np.shape(b)[1],np.shape(b)[2])    
    else :
        b = b.reshape(np.shape(b)[0],np.shape(b)[1])
    
    # Return the rotated data with correct 
    # dimensions.
    return b
