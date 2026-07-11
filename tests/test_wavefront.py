import pytest
import numpy as np
import lentil


def test_default_wavefront():
    w = lentil.Wavefront(wavelength=500e-9)
    assert np.array_equal(w.field, 1+0j)
    assert np.array_equal(w.intensity, 1)

def test_wavefront_rmul():
    w = lentil.Wavefront(wavelength=500e-9)
    p = lentil.Plane()
    assert p * w

def test_wavefront_propagate_image_non_pupil():
    w = lentil.Wavefront(wavelength=500e-9)
    with pytest.raises(TypeError):
        lentil.propagate_dft(w, shape=(64,64), pixelscale=5e-6)


@pytest.mark.parametrize('field_shape, field_shift, output_shape', [
    ((5,5), (-25,0), (10,10)),
    ((5,5), (25,0), (10,10)),
    ((5,5), (0,-25), (10,10)),
    ((5,5), (0,25), (10,10))
])
def test_overlap(field_shape, field_shift, output_shape):
    assert lentil.wavefront._overlap(field_shape, field_shift, output_shape) is False


def test_wavefront_default_state():
    w = lentil.Wavefront(wavelength=500e-9)
    assert w.z == 0
    assert w.pilot is None
    assert w.reference == 'planar'
    assert w.path == 0.0
    assert w.focal_length is None
    assert w.z_focus is None


def test_wavefront_focal_length_derived():
    w = lentil.Wavefront(wavelength=500e-9, focal_length=2, z=1)
    assert w.z_focus == 3
    assert w.focal_length == 2

    # z_focus is fixed; focal_length is derived from the current position
    w.z = 2
    assert w.z_focus == 3
    assert w.focal_length == 1


def test_wavefront_focal_length_setter():
    w = lentil.Wavefront(wavelength=500e-9, z=1)
    w.focal_length = 5
    assert w.z_focus == 6
    assert w.focal_length == 5

    w.focal_length = None
    assert w.z_focus is None
    assert w.focal_length is None

    w.focal_length = np.inf
    assert w.focal_length is None


def test_wavefront_pilot_from_diameter():
    w = lentil.Wavefront(wavelength=500e-9, diameter=1e-2, z=3)
    assert w.pilot is not None
    assert w.pilot.waist_radius == 5e-3
    assert w.pilot.waist_position == 3
    assert w.pilot.wavelength == 500e-9


def test_wavefront_pilot_explicit():
    from lentil.fresnel import GaussianBeam
    pilot = GaussianBeam(1e-3, 500e-9, waist_position=-1)
    w = lentil.Wavefront(wavelength=500e-9, diameter=1e-2, pilot=pilot)
    assert w.pilot is pilot


def test_wavefront_pilot_inherited_from_plane():
    amp = lentil.circle((64, 64), 24)
    p = lentil.Plane(amplitude=amp, pixelscale=1e-3)
    w = lentil.Wavefront(wavelength=500e-9)
    assert w.pilot is None

    w2 = w * p
    assert w2.pilot is not None
    assert np.isclose(w2.pilot.waist_radius, p.diameter/2)
    assert w2.pilot.waist_position == w.z


def test_wavefront_pilot_not_inherited_without_diameter():
    # a default Plane has no diameter to offer
    w = lentil.Wavefront(wavelength=500e-9) * lentil.Plane()
    assert w.pilot is None


def test_wavefront_state_preserved_through_multiply():
    amp = lentil.circle((64, 64), 24)
    p = lentil.Pupil(amplitude=amp, pixelscale=1e-3, focal_length=2)
    w = lentil.Wavefront(wavelength=500e-9, z=1)
    w2 = w * p
    assert w2.z == 1
    # a curvature-bearing plane produces a spherical-reference wavefront:
    # the (never sampled) ideal quadratic phase is absorbed by the
    # reference surface centered on the analytic focus
    assert w2.reference == 'spherical'
    assert w2.focal_length == 2
    assert w2.z_focus == 3
