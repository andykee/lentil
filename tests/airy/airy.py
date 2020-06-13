import numpy as np
from scipy import special


def airy(diameter, focal_length, wavelength, pixelscale, length, oversample=1):
    # https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation
    f_number = focal_length/diameter
    length *= oversample
    
    c = (length-1)/2
    x = np.arange(length, dtype=np.float)
    x -= c
    
    q = x * (pixelscale/oversample)
    X = (np.pi*q)/(wavelength*f_number)
   
    # if length is odd, the center value will be zero which will throw a 
    # divide by zero error. To avoid this, we'll set any zeros to machine
    # epsilon (eps)
    X[X == 0] = np.finfo(X.dtype).eps
    
    return (2*special.jn(1, X)/X)**2


def airy2(diameter, focal_length, wavelength, pixelscale, shape, oversample=1):
    # https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation
    
    f_number = focal_length/diameter
    
    shape = np.asarray(shape)
    shape *= oversample
    
    c = (shape-1)/2
    y, x = np.indices(shape, dtype=np.float)
    
    x -= c[1]
    y -= c[0]
    
    x *= (pixelscale/oversample)
    y *= (pixelscale/oversample)
    
    q = np.sqrt(x**2 + y**2)
    X = (np.pi*q)/(wavelength*f_number)
    
    # if length is odd, the center value will be zero which will throw a 
    # divide by zero error. To avoid this, we'll set any zeros to machine
    # epsilon (eps)
    X[X == 0] = np.finfo(X.dtype).eps
    
    return (2*special.jn(1, X)/X)**2
