import copy
from itertools import combinations
from operator import itemgetter

import numpy as np


def _standardize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def _standardize_bandpass(vec, default=()):
    if vec is None:
        vec = default
    vec = np.asarray(vec)
    if vec.shape == ():
        vec = vec[np.newaxis, ...]
    return vec


class _iterate_planes:
    def __init__(self, planes):
        self.planes = planes
        self.length = len(planes)
        self.n = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.n < self.length-1:
            plane = self.planes[self.n]
            next_plane = self.planes[self.n+1]
            self.n += 1
            return plane, next_plane
        else:
            raise StopIteration()


class _iterate_planes_reverse:
    def __init__(self, planes):
        self.planes = planes
        self.length = len(planes)
        self.n = 1

    def __iter__(self):
        return self

    def __next__(self):
        if self.n < self.length:
            plane = self.planes[-(self.n+1)]
            next_plane = self.planes[-self.n]
            self.n += 1
            return plane, next_plane
        else:
            raise StopIteration()


