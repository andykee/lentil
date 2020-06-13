import numpy as np

from lentil.wavefront import Wavefront
from lentil.util import centroid, circle
from tests.airy import airy2


diameter = 1
focal_length = 10
pupil_pixelscale = diameter/512
amplitude = circle((512, 512), 256)
wavelength = 650e-9

image_pixelscale = 5e-6
npix = 128
oversample = 10


def test_propagate_dft_airy():

    # we set oversample = 1 and apply an "oversampling" factor of 25 to pixelscale
    # so that we can more easily enforce an odd npix. having an odd npix enables
    # comparison between numerical propagations via dft2() and the airy2() function
    # since their DC pixels will be coincident.
    image_pixelscale = 5e-6/25
    npix = 511
    oversample = 1

    # numerical result
    w = Wavefront(wavelength, shape=amplitude.shape, pixelscale=pupil_pixelscale)
    w.data *= amplitude
    w.focal_length = focal_length
    w.propagate_dft(image_pixelscale, npix, oversample)
    psf = np.abs(w.data[0])**2
    psf = psf/np.max(psf)

    # analytical result
    airy = airy2(diameter, focal_length, wavelength, image_pixelscale, (npix, npix), oversample)
    airy = airy/np.max(airy)

    assert np.all(np.isclose(psf, airy, atol=1e-3))


def test_propagate_dft_shift():
    shift = np.random.uniform(low=-25, high=25, size=2)
    w = Wavefront(wavelength, shape=amplitude.shape, pixelscale=pupil_pixelscale)
    w.data *= amplitude
    w.focal_length = focal_length
    w.propagate_dft(image_pixelscale, npix, oversample, shift)

    psf = np.abs(w.data[0])**2
    psf = psf/np.max(psf)

    psf[psf < 1e-2] = 0

    center = np.array([(npix*oversample)//2, (npix*oversample)//2])
    psf_centroid = np.asarray(centroid(psf))
    psf_shift = psf_centroid - center

    assert np.all(np.isclose(shift, psf_shift, atol=1e-2))

# TODO: get rid of these unused (and covered in test_propagate) tests
# def test_propagate_dft_int_shift():
#     shift = np.random.randint(low=-100, high=100, size=2)
#     shift_m = shift * (image_pixelscale/oversample)
#
#     w = Wavefront(wavelength, shape=amplitude.shape, pixelscale=pupil_pixelscale)
#     w.data *= amplitude
#     w.focal_length = focal_length
#     w.tilt = [Shift(x=shift_m[0], y=shift_m[1])]
#     w.propagate_dft(image_pixelscale, npix, oversample)
#
#     psf = np.abs(w.data[0])**2
#     psf = psf/np.max(psf)
#
#     center = np.array([(npix*oversample)//2, (npix*oversample)//2])
#     psf_centroid = np.asarray(centroid(psf))
#     psf_shift = psf_centroid - center
#
#     assert np.all(psf_shift < 1e-8)
#     assert np.all(np.equal(shift, w.shift))
#
#
# def test_propagate_dft_float_shift():
#     shift = np.random.uniform(low=-100, high=100, size=2)
#     shift_m = shift * (image_pixelscale/oversample)
#
#     w = Wavefront(wavelength, shape=amplitude.shape, pixelscale=pupil_pixelscale)
#     w.data *= amplitude
#     w.focal_length = focal_length
#     w.tilt = [Shift(x=shift_m[0], y=shift_m[1])]
#     w.propagate_dft(image_pixelscale, npix, oversample)
#
#     psf = np.abs(w.data[0])**2
#     psf = psf/np.max(psf)
#
#     fix_shift = np.fix(shift)
#     res_shift = shift - fix_shift
#
#     center = np.array([(npix*oversample)//2, (npix*oversample)//2])
#     psf_shifted_center = center + res_shift
#     psf_centroid = np.asarray(centroid(psf))
#     psf_shift_error = psf_centroid - psf_shifted_center
#
#     assert np.all(psf_shift_error < 1e-3)
#     assert np.all(np.equal(fix_shift, w.shift))
#
#
# def test_propagate_dft_angle():
#     shift = np.random.uniform(low=-100, high=100, size=2)
#     angle = shift * (image_pixelscale/oversample) / focal_length
#
#     w = Wavefront(wavelength, shape=amplitude.shape, pixelscale=pupil_pixelscale)
#     w.data *= amplitude
#     w.focal_length = focal_length
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft(image_pixelscale, npix, oversample)
#
#     psf = np.abs(w.data[0])**2
#     psf = psf/np.max(psf)
#
#     fix_shift = np.fix(shift)
#     res_shift = shift - fix_shift
#
#     center = np.array([(npix*oversample)//2, (npix*oversample)//2])
#     psf_shifted_center = center + res_shift
#     psf_centroid = np.asarray(centroid(psf))
#     psf_shift_error = psf_centroid - psf_shifted_center
#
#     assert np.all(psf_shift_error < 1e-3)
#     assert np.all(np.equal(fix_shift, w.shift))


# def test_propagate_dft_shift():
#     shift = np.random.uniform(low=-100, high=100, size=2)
#     angle = (shift*pupil_pixelscale)/focal_length
#
#     # propagate_dft result
#     w = Wavefront(wavelength, shape=amplitude.shape, pixelscale=pupil_pixelscale)
#     w.data *= amplitude
#     w.focal_length = focal_length
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft(image_pixelscale, npix, oversample)
#     psf = np.abs(w.data[0])**2
#     psf = psf/np.max(psf)
#
#     # threshold the PSF so that the centroid calculation makes sense
#     psf[psf < 0.5] = 0
#
#     center = np.array([(npix*oversample)//2, (npix*oversample)//2])
#     psf_centroid = np.asarray(centroid(psf))
#     psf_shift = psf_centroid - center
#     assert np.all(np.isclose(shift, psf_shift, atol=0.1))
#
#
# def test_propagate_dft_coarse_int_shift():
#     shift = np.round(np.random.uniform(low=-100, high=100, size=2))
#     angle = (shift*pixelscale)/p.focal_length
#
#     # propagate_dft_coarse result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft_coarse(pixelscale, npix, npix, oversample)
#     psf = np.abs(w.data)**2
#     psf = psf/np.max(psf)
#
#     # threshold the PSF so that the centroid calculation makes sense
#     psf[psf < 0.5] = 0
#
#     center = np.array([(npix*oversample)//2, (npix*oversample)//2])
#     psf_centroid = np.asarray(centroid(psf))
#     psf_shift = psf_centroid - center
#     assert np.all(np.isclose(shift, psf_shift, atol=0.1))
#
#
# def test_propagate_dft_coarse_float_shift():
#     shift = np.random.uniform(low=-100, high=100, size=2)
#     angle = (shift*pixelscale)/p.focal_length
#
#     # propagate_dft_coarse result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft_coarse(pixelscale, npix, npix, oversample)
#     psf = np.abs(w.data)**2
#     psf = psf/np.max(psf)
#
#     # threshold the PSF so that the centroid calculation makes sense
#     psf[psf < 0.5] = 0
#
#     center = np.array([(npix*oversample)//2, (npix*oversample)//2])
#     psf_centroid = np.asarray(centroid(psf))
#     psf_shift = psf_centroid - center
#     assert np.all(np.isclose(shift, psf_shift, atol=0.1))
#
#
# def test_propagate_dft_coarse_huge_shift():
#     # A shift so large the chip falls entirely out of the desired area
#     shift = np.asarray([20000, 20000])
#     angle = (shift*pixelscale)/p.focal_length
#
#     psf = np.zeros([npix, npix], dtype=np.complex128)
#
#     # propagate_dft_coarse result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft_coarse(pixelscale, npix, npix, oversample)
#     psf_coarse = np.abs(w.data)**2
#     # This should return zeros so top avoid a divide by zero warning, we won't normalize
#
#     assert np.all(np.isclose(psf, psf_coarse, atol=5e-4))
#
#
# def test_propagate_dft_vs_propagate_dft_coarse_no_shift():
#
#     # propagate_dft result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.propagate_dft(pixelscale, npix, oversample)
#     psf = np.abs(w.data)**2
#     psf = psf/np.max(psf)
#
#     # propagate_dft_coarse result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.propagate_dft_coarse(pixelscale, npix, npix, oversample)
#     psf_coarse = np.abs(w.data)**2
#     psf_coarse = psf_coarse/np.max(psf_coarse)
#
#     assert np.all(np.isclose(psf, psf_coarse, atol=5e-4))
#
#
# def test_propagate_dft_vs_propagate_dft_coarse_shift():
#     shift = np.random.uniform(low=-100, high=100, size=2)
#     angle = (shift*pixelscale)/p.focal_length
#
#     # propagate_dft result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft(pixelscale, npix, oversample)
#     psf = np.abs(w.data)**2
#     psf = psf/np.max(psf)
#
#     # propagate_dft_coarse result
#     w = Wavefront(wavelength, shape=p.shape, pixelscale=p.pixelscale)
#     w = p.multiply(w)
#     w.tilt = [Angle(x=angle[0], y=angle[1])]
#     w.propagate_dft_coarse(pixelscale, npix, npix, oversample)
#     psf_coarse = np.abs(w.data)**2
#     psf_coarse = psf_coarse/np.max(psf_coarse)
#
#     assert np.all(np.isclose(psf, psf_coarse, atol=5e-4))
