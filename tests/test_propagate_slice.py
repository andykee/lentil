# These tests deal specifically with ensuring sliced and offset propagations
# are equivalent to their full-sized counterparts

import numpy as np
import lentil
import lentil.helper


def test_propagate_slice_one():
    # Single monolithic aperture

    diameter = 1
    focal_length = 20
    wavelength = 500e-9
    du = 5e-6
    dx = diameter/128

    shift = np.random.uniform(low=-50, high=50, size=2).astype(int)
    amp = lentil.circle((256,256), 256//4, shift=shift)

    rho, theta = lentil.zernike_coordinates(amp, shift=shift)
    coeffs = np.random.uniform(low=-1, high=1, size=11)*100e-9
    coeffs[0:3] = 0 # no ptt
    opd = lentil.zernike_compose(amp, coeffs, rho=rho, theta=theta)

    alpha = (dx*du)/(wavelength*focal_length)
    oversample = 5
    shape = 64

    phasor = amp*np.exp(-2*1j*np.pi*opd/wavelength)
    F = lentil.fourier.dft2(phasor, alpha=alpha/oversample, shape=shape*oversample)


    slc = lentil.helper.boundary_slice(amp)
    ofst = lentil.helper.slice_offset(slc, shape=amp.shape)

    phasor = amp[slc]*np.exp(-2*1j*np.pi*opd[slc]/wavelength)
    F_slc = lentil.fourier.dft2(phasor, alpha=alpha/oversample, shape=shape*oversample, offset=ofst)

    assert np.allclose(F, F_slc)


def test_propagate_slice_multi():
    # Segmented aperture
    diameter = 1
    focal_length = 20
    wavelength = 500e-9
    du = 5e-6
    dx = diameter/102

    n = 256

    amp1 = lentil.circle((n,n), n//5, shift=(0, -0.3*n))
    slc1 = lentil.helper.boundary_slice(amp1)
    ofst1 = lentil.helper.slice_offset(slc1, shape=amp1.shape)

    rho, theta = lentil.zernike_coordinates(amp1, shift=(0, -0.3*n))
    coeffs = np.random.uniform(low=-1, high=1, size=11)*100e-9
    coeffs[0:3] = 0 # no ptt
    opd1 = lentil.zernike_compose(amp1, coeffs, rho=rho, theta=theta)


    amp2 = lentil.circle((n,n), n//5, shift=(0, .3*n))
    slc2 = lentil.helper.boundary_slice(amp2)
    ofst2 = lentil.helper.slice_offset(slc2, shape=amp2.shape)

    rho, theta = lentil.zernike_coordinates(amp2, shift=(0, 0.3*n))
    coeffs = np.random.uniform(low=-1, high=1, size=11)*100e-9
    coeffs[0:3] = 0 # no ptt
    opd2 = lentil.zernike_compose(amp2, coeffs, rho=rho, theta=theta)

    alpha = (dx*du)/(wavelength*focal_length)
    oversample = 5
    shape = 64

    amp = amp1 + amp2
    opd = opd1 + opd2

    phasor = amp*np.exp(-2*1j*np.pi*opd/wavelength)
    F = lentil.fourier.dft2(phasor, alpha=alpha/oversample, shape=shape*oversample)

    phasor1 = amp1[slc1]*np.exp(-2*1j*np.pi*opd1[slc1]/wavelength)
    F1 = lentil.fourier.dft2(phasor1, alpha=alpha/oversample, shape=shape*oversample, offset=ofst1)

    phasor2 = amp2[slc2]*np.exp(-2*1j*np.pi*opd2[slc2]/wavelength)
    F2 = lentil.fourier.dft2(phasor2, alpha=alpha/oversample, shape=shape*oversample, offset=ofst2)

    F_slc = F1 + F2

    np.allclose(F, F_slc)
