import numpy as np
from scipy import interpolate

from tests.airy import airy, airy2


def test_airy():
    diameter = 1
    focal_length = 10
    f_number = focal_length/diameter
    wave = 650e-9
    pixelscale = 5e-6/1000  # divide by 1000 here to fake oversampling
    length = 2**14  # an appropriately large number
    oversample = 1

    airy1 = airy(diameter, focal_length, wave, pixelscale, length, oversample)

    c = (length-1)/2
    x = np.arange(length, dtype=np.float)
    x -= c
    x *= pixelscale

    airy_interp = interpolate.interp1d(x, airy1)

    # TODO: move this in to the documentation and just point to that as a reference here
    # derivation of airy FWHM expression 
    # ----------------------------------
    # the half-maximum of the central airy disk occurs at x_hwhm ~= 1.61633 [1]
    # since x = (pi*q)/(wave*f_number), q = (x * wave * f_number)/pi where
    # q is the radial distance from the optics axis in the observation
    # (or focal) plane
    #
    # substituting in x_hwhm,
    # q_hwhm = (1.61633/pi) * wave * f_number
    #
    # because we are more interested in FWHM, we multiply by 2
    # q_fwhm = 2*(1.61633/pi) * wave * f_number = 1.028987 * wave * f_number
    #
    # Refs:
    # [1] https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation

    # verify the FWHM occurs at 1.029 * wave * f_number
    j1_hwhm = 1.61633
    fwhm = 2*(j1_hwhm/np.pi)*wave*f_number  # actual linear distance
    assert(np.abs(airy_interp(fwhm/2) - 0.5) < 1e-5)

    # verify the first five zeros (as listed on Wikipedia) occur in the
    # right places
    j1_zeros = np.asarray([3.8317, 7.0156, 10.1735, 13.3237, 16.4706])
    for j1_zero in j1_zeros:
        assert(airy_interp(j1_zero/np.pi * wave * f_number) < 1e-7)


def test_airy2():
    # evaluate the 2D Airy code against the 1D Airy code. since we've done 
    # some pretty exhaustive testing of the 1D Airy code, this simple 
    # comparison should be sufficient
    
    diameter = 1
    focal_length = 10
    wave = 650e-9
    pixelscale = 5e-6/10  # divide by 10 here to fake oversampling
    npix = 511
    oversample = 1
    slice_index = npix//2
    
    a1 = airy(diameter, focal_length, wave, pixelscale, npix, oversample)
    a2 = airy2(diameter, focal_length, wave, pixelscale, (npix, npix), oversample)

    res = a2[:, slice_index] - a1
    assert(np.all(np.abs(res) < 1e-9))
