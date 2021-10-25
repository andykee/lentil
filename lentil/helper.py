# Lentil internal helper functions

import numpy as np


def sanitize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def sanitize_bandpass(vec, default=()):
    if vec is None:
        vec = default
    vec = np.asarray(vec)
    if vec.shape == ():
        vec = vec[np.newaxis, ...]
    return vec


def dft_alpha(dx, du, wave, z, oversample):
    return (dx*du)/(wave*z*oversample)
