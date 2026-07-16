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
    assert p.pixelscale is None
    assert np.all(p.amplitude == 1)
    assert np.all(p.opd == 0)
    assert p.mask == p.amplitude


def test_amp_alias():
    p = lentil.Plane(amp=10)
    assert p.amplitude == 10


def test_amp_alias_error():
    with pytest.raises(AttributeError):
        lentil.Plane(amplitude=10, amp=10)


def test_plane_overload_methods():
    class Plane(lentil.Plane):
        def __init__(self):
            super().__init__()

        def __amp__(self):
            return np.array(1)

        def __opd__(self):
            return np.array(2)

        def __mask__(self):
            return np.array(3)

    p = Plane()

    assert p.amplitude == 1
    assert p.opd == 2
    assert p.mask == 3


def test_plane_opd_overload_call_mask():
    class Plane(lentil.Plane):

        def __opd__(self):
            return 2 * self.mask

        def __mask__(self):
            return 3

    p = Plane()
    assert p.opd == 6


def test_plane_fit_tilt_inplace():
    p = RandomPlane()
    p_copy = p.fit_tilt(inplace=False)
    p_inplace = p.fit_tilt(inplace=True)

    assert p_copy is not p
    assert p_inplace is p


class ComputedOPDPupil(lentil.Pupil):
    # opd is computed on every access via __opd__()
    def __init__(self, coeffs, **kwargs):
        super().__init__(**kwargs)
        self.coeffs = np.asarray(coeffs)

    def __opd__(self):
        return lentil.zernike_compose(self.mask, self.coeffs)


def _fit_residual_tilt(plane):
    ptt = plane.ptt_vector
    t = np.linalg.lstsq(ptt.T, np.asarray(plane.opd).ravel(), rcond=None)[0]
    return t[1], t[2]


def test_fit_tilt_computed_opd_removes_tilt():
    mask = lentil.circle((64, 64), 30)
    coeffs = [0, 300e-9, 200e-9, 100e-9]  # piston, x-tilt, y-tilt, focus
    p = ComputedOPDPupil(coeffs, amplitude=mask, mask=mask, pixelscale=1/64,
                         focal_length=10)

    pf = p.fit_tilt(inplace=False)

    # tilt is actually removed from the returned opd
    tx, ty = _fit_residual_tilt(pf)
    assert np.isclose(tx, 0, atol=1e-12)
    assert np.isclose(ty, 0, atol=1e-12)

    # bookkept exactly once
    assert len(pf.tilt) == 1

    # the result is a frozen snapshot
    assert pf.frozen

    # and the original is untouched (fit_tilt returned a copy)
    assert not p.frozen
    otx, oty = _fit_residual_tilt(p)
    assert not np.isclose(otx, 0, atol=1e-12)


def test_fit_tilt_idempotent_on_frozen():
    mask = lentil.circle((64, 64), 30)
    coeffs = [0, 300e-9, 200e-9, 100e-9]
    p = ComputedOPDPupil(coeffs, amplitude=mask, mask=mask, pixelscale=1/64,
                         focal_length=10).fit_tilt(inplace=False)

    # already frozen -> fit_tilt must not raise and must leave a tilt-free opd
    pf = p.fit_tilt(inplace=False)
    tx, ty = _fit_residual_tilt(pf)
    assert np.isclose(tx, 0, atol=1e-12)
    assert np.isclose(ty, 0, atol=1e-12)


def test_wavefront_plane_mul():
    p = RandomPlane()
    w = lentil.Wavefront(650e-9)

    w1 = w * p

    slc = lentil.helper.boundary_slice(p.mask)
    phasor = p.amplitude[slc] * np.exp(2*np.pi*1j*p.opd[slc]/w.wavelength)

    assert np.array_equal(w1.data[0].data, phasor)


def test_wavefront_plane_rmul():
    p = RandomPlane()
    w = lentil.Wavefront(650e-9)

    w1 = p * w

    slc = lentil.helper.boundary_slice(p.mask)
    phasor = p.amplitude[slc] * np.exp(2*np.pi*1j*p.opd[slc]/w.wavelength)

    assert np.array_equal(w1.data[0].data, phasor)


def test_wavefront_plane_imul():
    p = RandomPlane()
    w = lentil.Wavefront(650e-9)

    w *= p

    slc = lentil.helper.boundary_slice(p.mask)
    phasor = p.amplitude[slc] * np.exp(2*np.pi*1j*p.opd[slc]/w.wavelength)

    assert np.array_equal(w.data[0].data, phasor)


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
    w = p * w
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
    x, y = dt.__shift__(wavelength=650e-9, xs=x0, ys=y0)
    assert np.all((x == x0, y == y0))


def test_dispersive_tilt_shift():
    dt = lentil.DispersiveTilt(trace=[1, 1], dispersion=[1, 650e-9])
    wave = 900e-9
    x, y = dt.__shift__(wavelength=wave)

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


@pytest.mark.parametrize('plane, ptype', [
    (lentil.Plane, lentil.ptype(None)),
    (lentil.Pupil, lentil.pupil),
    (lentil.Image, lentil.image),
    (lentil.plane._TiltBase, lentil.tilt)
])
def test_default_ptype(plane, ptype):
    p = plane()
    assert p.ptype == ptype


def test_tilt_ptype_overload():
    t = lentil.plane._TiltBase(ptype=lentil.pupil)
    assert t.ptype == lentil.pupil


def test_dispersivetilt_overload():
    class Plane(lentil.DispersiveTilt):
        dispersion = [1, 0]
        trace = [1, 0]

    p = Plane()
    assert np.allclose(p.__shift__(1), np.sqrt(2)/2)
