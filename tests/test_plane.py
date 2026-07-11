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


def test_mul_ptype_none_preserves_conjugate():
    # a ptype-none plane (screen, mask, lens) does not move the wavefront
    # off a conjugate: pupil x none -> pupil, image x none -> image
    amp = lentil.circle((32, 32), 12)
    screen = lentil.Plane(amplitude=amp, pixelscale=1e-3)

    w = lentil.Wavefront(wavelength=500e-9, ptype=lentil.pupil)
    assert (w * screen).ptype == lentil.pupil

    w = lentil.Wavefront(wavelength=500e-9, ptype=lentil.image)
    assert (w * screen).ptype == lentil.image

    w = lentil.Wavefront(wavelength=500e-9)
    assert (w * screen).ptype == lentil.none


def test_lens_plane_wave():
    w = lentil.Wavefront(wavelength=500e-9)
    w2 = w * lentil.Lens(focal_length=2)
    assert w2.focal_length == 2
    assert w2.z_focus == 2
    assert w2.ptype == lentil.none


def test_lens_law_combination():
    # 1/fl_out = 1/fl_in + 1/f
    w = lentil.Wavefront(wavelength=500e-9, focal_length=2)
    w2 = w * lentil.Lens(focal_length=2)
    assert np.isclose(w2.focal_length, 1)


def test_lens_collimates():
    # a wavefront diverging from a virtual focus one focal length behind
    # the lens is collimated by it
    w = lentil.Wavefront(wavelength=500e-9, focal_length=-2)
    w2 = w * lentil.Lens(focal_length=2)
    assert w2.focal_length is None


def test_lens_no_power():
    w = lentil.Wavefront(wavelength=500e-9, focal_length=3)
    w2 = w * lentil.Lens(focal_length=np.inf)
    assert w2.focal_length == 3


def test_lens_stack_combined_power():
    # two thin lenses in contact combine their powers
    w = lentil.Wavefront(wavelength=500e-9)
    w2 = w * lentil.Lens(focal_length=1) * lentil.Lens(focal_length=1)
    assert np.isclose(w2.focal_length, 0.5)


def test_lens_updates_pilot():
    # collimated pilot focused by a lens comes to a waist ~one focal
    # length away with the diffraction-limited radius
    wavelength, diameter, f = 1e-6, 0.01, 1
    w = lentil.Wavefront(wavelength=wavelength, diameter=diameter)
    w2 = w * lentil.Lens(focal_length=f)
    assert np.isclose(w2.pilot.waist_position, f, rtol=1e-3)
    assert np.isclose(w2.pilot.waist_radius,
                      wavelength*f/(np.pi*diameter/2), rtol=1e-3)
    # the input wavefront's pilot is unchanged
    assert w.pilot.waist_position == 0


def test_lens_never_samples_quadratic_phase():
    # only aberration OPD is sampled into the field; the ideal quadratic
    # phase lives in the analytic curvature state
    amp = lentil.circle((32, 32), 12)
    opd = 1e-7 * lentil.circle((32, 32), 12)
    lens = lentil.Lens(focal_length=2, amplitude=amp, opd=opd,
                       pixelscale=1e-3)
    w = lentil.Wavefront(wavelength=500e-9)
    w2 = w * lens

    expected = amp * np.exp(2*np.pi*1j*opd/500e-9)
    assert np.allclose(w2.field, expected)
    assert w2.focal_length == 2


def test_lens_pupil_equivalence():
    # for a plane-wave input, Lens and Pupil produce the same curvature
    # state; they differ only in ptype
    amp = lentil.circle((32, 32), 12)
    lens = lentil.Lens(focal_length=2, amplitude=amp, pixelscale=1e-3)
    pupil = lentil.Pupil(focal_length=2, amplitude=amp, pixelscale=1e-3)

    wl = lentil.Wavefront(wavelength=500e-9) * lens
    wp = lentil.Wavefront(wavelength=500e-9) * pupil

    assert wl.focal_length == wp.focal_length
    assert wl.z_focus == wp.z_focus
    assert np.isclose(wl.pilot.waist_position, wp.pilot.waist_position)
    assert np.isclose(wl.pilot.waist_radius, wp.pilot.waist_radius)
    assert wl.ptype == lentil.none
    assert wp.ptype == lentil.pupil


def test_pupil_updates_pilot():
    amp = lentil.circle((32, 32), 12)
    pupil = lentil.Pupil(focal_length=2, amplitude=amp, pixelscale=1e-3)
    w = lentil.Wavefront(wavelength=500e-9) * pupil
    assert w.pilot is not None
    assert np.isclose(w.pilot.waist_position, 2, rtol=1e-3)


def test_lens_zero_focal_length():
    with pytest.raises(ValueError):
        lentil.Lens(focal_length=0)
