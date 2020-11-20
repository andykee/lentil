import numpy as np


def expc(x):
    """Computes np.exp(1.0j * x) for real valued x... FAST!"""
    x = np.asarray(x)
    if x.dtype == np.float32:
        out = np.empty(x.shape, dtype=np.complex64)
    elif x.dtype == np.float64:
        out = np.empty(x.shape, dtype=np.complex128)
        
    out.real = np.cos(x)
    out.imag = np.sin(x)
    return out
