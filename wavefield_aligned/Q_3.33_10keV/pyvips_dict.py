import numpy as np
import pyvips

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

def pyvips_rotate(x,angle_):       
    angle_rad = angle_*(np.pi/180)
    tc = np.cos(angle_rad)
    ts = np.sin(angle_rad)

    rows, cols = x.shape[0], x.shape[1]
    c_ = np.array((rows,cols))/2 - 0.5
    
    if len(x.shape)!=3 : 
        height, width = x.shape
        bands = 1
    else :
        height, width,bands = x.shape
    
    im = pyvips.Image.new_from_memory(x.reshape(width * height * bands).data, width, height, bands,
                                      dtype_to_format[str(x.dtype)])

    inter = pyvips.vinterpolate.Interpolate.new('bilinear')

    im = im.affine([tc, ts, -ts, tc], interpolate = inter,
                   idx=-c_[1],idy=-c_[0],
                   odx=c_[1],ody=c_[0],
                   oarea=[0, 0, im.width, im.height])

    b = np.ndarray(buffer=im.write_to_memory(),dtype=format_to_dtype[im.format],shape=[im.height, im.width, im.bands])
    del im
    if len(x.shape)==3 : 
        b = b.reshape(np.shape(b)[0],np.shape(b)[1],np.shape(b)[2])    
    else :
        b = b.reshape(np.shape(b)[0],np.shape(b)[1])
    
    return b