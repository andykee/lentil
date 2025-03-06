import math
import numpy as np

import pytest
import lentil


size = (10, 10)


class RandomPlane(lentil.Plane):
    def __init__(self):
        super().__init__(pixelscale=1,
                         amplitude=np.random.uniform(size=size),
                         opd=np.random.uniform(size=size))


def test_default_plane():
    # Ensure that a default Plane creates an object that won't have any
    # impact on an optical system (a perfect optic with no wavefront error 
    # and perfect optical and spectral transmission).

    p = lentil.Plane()
    assert p.pixelscale == None
    assert np.all(p.amplitude == 1)
    assert np.all(p.opd == 0)
    assert p.mask == p.amplitude


def test_amp_alias():
    p = lentil.Plane(amp=10)
    assert p.amplitude == 10


def test_amp_alias_error():
    with pytest.raises(TypeError):
        lentil.Plane(amplitude=10, amp=10)


class PlaneAttributeOverload(lentil.Plane):
    def __init__(self):
        super().__init__()

    @property
    def amplitude(self):
        return np.array(1)
    
    @property
    def opd(self):
        return np.array(2)
    

def test_plane_overload_propertyes():
    p = PlaneAttributeOverload()
    
    assert p.amplitude == 1
    assert p.opd == 2
    assert p.mask == 1

    with pytest.raises(AttributeError):
        p.opd = 5

    with pytest.raises(AttributeError):
        p.amplitude = 5


def test_plane_fit_tilt_inplace():
    p = RandomPlane()
    p_copy = p.fit_tilt(inplace=False)
    p_inplace = p.fit_tilt(inplace=True)

    assert p_copy is not p
    assert p_inplace is p


def test_wavefront_plane_multiply():
    p = RandomPlane()
    w = lentil.Wavefront(650e-9)

    w1 = p.multiply(w)

    slc = lentil.helper.boundary_slice(p.mask)
    phasor = p.amplitude[slc] * np.exp(2*np.pi*1j*p.opd[slc]/w.wavelength)

    assert np.array_equal(w1.data[0].data, phasor)


def test_wavefront_plane_multiply_overlapping_segment_slices():
    seg = lentil.hexagon((64, 64), 32, antialias=False)
    seg = seg[5:60, :]

    segmask = np.zeros((2, 128, 128))
    segmask[0, 0:55, 2:66] = seg
    segmask[1, 29:84, 55:119] = seg
    mask = np.sum(segmask, axis=0)

    pupil = lentil.Pupil(amplitude=mask, mask=segmask, pixelscale=1 / 256, focal_length=10)

    w = lentil.Wavefront(500e-9)
    w *= pupil

    assert np.array_equal(mask, w.intensity)


class CircularPupil(lentil.Pupil):
    def __init__(self):
        super().__init__(focal_length=10,
                         pixelscale=2/256,
                         amplitude=lentil.circle((256, 256), 128),
                         opd=np.zeros((256, 256)))


def test_wavefront_pupil_multiply():
    p = CircularPupil()
    w = lentil.wavefront.Wavefront(650e-9)
    w = p.multiply(w)
    phasor = p.amplitude * np.exp(1j*p.opd * 2 * np.pi / w.wavelength)

    assert np.array_equal(w.data[0].data, phasor)
    assert w.focal_length == p.focal_length


def test_pupil_rescale_power():
    p = CircularPupil()
    pr = p.rescale(3)

    amp_power = np.sum(np.abs(p.amplitude)**2)
    ampr_power = np.sum(np.abs(pr.amplitude)**2)
    assert math.isclose(amp_power, ampr_power, rel_tol=1e-3)


def test_dispersive_tilt_center():
    dispersion = [1, 650e-9]
    trace = [2, 0]
    dt = lentil.DispersiveTilt(dispersion=dispersion, trace=trace)

    x0 = np.random.uniform(low=-5, high=5)
    y0 = np.random.uniform(low=-5, high=5)
    x, y = dt.shift(wavelength=650e-9, xs=x0, ys=y0)
    assert np.all((x == x0, y == y0))


def test_dispersive_tilt_shift():
    dt = lentil.DispersiveTilt(trace=[1,1], dispersion=[1,650e-9])
    wave = 900e-9
    x, y = dt.shift(wavelength=wave)

    assert x == (wave - dt.dispersion[1])/np.sqrt(2)
    assert y == 1+x


class FreezePlane(lentil.Plane):

    def __opd__(self):
        return np.random.uniform()


def test_freeze():
    p = FreezePlane()
    p.freeze()
    a = p.opd
    b = p.opd
    assert a == b
