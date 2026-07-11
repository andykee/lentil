import pytest
import numpy as np

from lentil.fresnel import GaussianBeam


def test_gaussian_beam_rayleigh_distance():
    # w0 = 1 mm, wavelength = 1 um -> z_R = pi meters
    b = GaussianBeam(1e-3, 1e-6)
    assert np.isclose(b.rayleigh_distance, np.pi)


def test_gaussian_beam_radius():
    b = GaussianBeam(1e-3, 1e-6, waist_position=2)
    assert b.radius(2) == 1e-3
    # beam expands by sqrt(2) at the Rayleigh distance
    assert np.isclose(b.radius(2 + b.rayleigh_distance), np.sqrt(2)*1e-3)
    assert np.isclose(b.radius(2 - b.rayleigh_distance), np.sqrt(2)*1e-3)


def test_gaussian_beam_phase_radius():
    b = GaussianBeam(1e-3, 1e-6, waist_position=2)
    zr = b.rayleigh_distance

    # infinite at the waist
    assert np.isinf(b.phase_radius(2))

    # maximum curvature (minimum |R|) of 2*z_R at the Rayleigh distance
    assert np.isclose(b.phase_radius(2 + zr), 2*zr)
    assert np.isclose(b.phase_radius(2 - zr), -2*zr)

    # diverging (R > 0) past the waist, converging (R < 0) before it
    assert b.phase_radius(2 + 10) > 0
    assert b.phase_radius(2 - 10) < 0


def test_gaussian_beam_gouy():
    b = GaussianBeam(1e-3, 1e-6)
    zr = b.rayleigh_distance
    assert np.isclose(b.gouy(zr), np.pi/4)
    assert np.isclose(b.gouy(-zr), -np.pi/4)
    assert b.gouy(0) == 0


def test_gaussian_beam_inside():
    b = GaussianBeam(1e-3, 1e-6, waist_position=1)
    zr = b.rayleigh_distance
    assert b.inside(1)
    assert b.inside(1 + zr)
    assert b.inside(1 - zr)
    assert not b.inside(1 + 1.01*zr)
    assert not b.inside(1 - 1.01*zr)


@pytest.mark.parametrize('z', [7.5, -3])
def test_gaussian_beam_from_radius_roundtrip(z):
    # waist finding inverts radius() and phase_radius() on both sides
    # of the waist
    b = GaussianBeam(1e-3, 1e-6, waist_position=2)
    c = GaussianBeam.from_radius(b.radius(z), b.phase_radius(z),
                                 b.wavelength, z=z)
    assert np.isclose(c.waist_radius, b.waist_radius)
    assert np.isclose(c.waist_position, b.waist_position)


def test_gaussian_beam_from_radius_at_waist():
    b = GaussianBeam.from_radius(1e-3, np.inf, 1e-6, z=4)
    assert b.waist_radius == 1e-3
    assert b.waist_position == 4


def test_gaussian_beam_lens_collimated():
    # a large collimated beam focused by a lens comes to a waist one
    # focal length away (geometric limit) with the diffraction-limited
    # waist radius w0' = lambda*f/(pi*w)
    w, wavelength, f = 0.01, 1e-6, 1
    b = GaussianBeam(w, wavelength, waist_position=0)
    c = b.lens(f, z=0)
    assert np.isclose(c.waist_position, f, rtol=1e-4)
    assert np.isclose(c.waist_radius, wavelength*f/(np.pi*w), rtol=1e-4)


def test_gaussian_beam_lens_imaging():
    # 2f-2f imaging: an object waist 2f before the lens is reimaged 2f
    # after the lens at unit magnification (in the z_R << f limit)
    wavelength, f = 1e-6, 0.1
    b = GaussianBeam(1e-6, wavelength, waist_position=0)
    c = b.lens(f, z=2*f)
    assert np.isclose(c.waist_position, 4*f, rtol=1e-3)
    assert np.isclose(c.waist_radius, b.waist_radius, rtol=1e-3)


def test_gaussian_beam_lens_collimating():
    # a lens with f = R placed one focal length past the waist
    # collimates the beam: new waist at the lens with w0' = w(z)
    b = GaussianBeam(1e-6, 1e-6, waist_position=0)
    z = 0.5
    c = b.lens(b.phase_radius(z), z=z)
    assert c.waist_position == z
    assert np.isclose(c.waist_radius, b.radius(z))


def test_gaussian_beam_lens_flat():
    # a lens with infinite focal length has no effect
    b = GaussianBeam(1e-3, 1e-6, waist_position=3)
    c = b.lens(np.inf, z=10)
    assert c.waist_radius == b.waist_radius
    assert c.waist_position == b.waist_position


def test_gaussian_beam_invalid():
    with pytest.raises(ValueError):
        GaussianBeam(-1, 1e-6)
    with pytest.raises(ValueError):
        GaussianBeam(1e-3, 0)
