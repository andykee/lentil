import numpy as np
import scipy

from lentil.util import *


def test_rebin():
    img = np.array([[1, 1, 2, 2], [1, 1, 2, 2], [3, 3, 4, 4], [3, 3, 4, 4]])
    factor = 2
    img_rebinned = rebin(img, factor)
    assert np.array_equal(img_rebinned, np.array([[4, 8], [12, 16]]))


def test_rebin_cube():
    img = np.zeros((2, 4, 4))
    img[0] = np.array([[1, 1, 2, 2], [1, 1, 2, 2], [3, 3, 4, 4], [3, 3, 4, 4]])
    img[1] = np.array([[2, 2, 4, 4], [2, 2, 4, 4], [6, 6, 8, 8], [6, 6, 8, 8]])
    factor = 2
    img_rebinned = rebin(img, factor)
    img_rebinned_expected = np.zeros((2, 2, 2))
    img_rebinned_expected[0] = np.array([[4, 8], [12, 16]])
    img_rebinned_expected[1] = np.array([[8, 16], [24, 32]])
    assert np.array_equal(img_rebinned, img_rebinned_expected)


def test_mesh_nonsquare():
    xx, yy = mesh((256, 512))
    assert xx.shape == (256, 512)
