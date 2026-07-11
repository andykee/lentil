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


import lentil


def _gaussian_wavefront(wl=1e-6, dx=4e-6, n=512, w0_px=40):
    w0 = w0_px*dx
    f = np.fft.fftshift(np.fft.fftfreq(n, 1/(n*dx)))
    r, c = np.meshgrid(f, f, indexing='ij')
    E0 = np.exp(-(r**2 + c**2)/w0**2)
    w = lentil.Wavefront(wl, pixelscale=dx)
    w = w * lentil.Plane(amplitude=E0, pixelscale=dx)
    return w, GaussianBeam(w0, wl)


def test_ptp_gaussian_closed_form():
    # a Gaussian propagated with the angular spectrum kernel matches the
    # analytic Gaussian beam solution:
    #   E(r, z) = (w0/w) exp(-r^2/w^2) exp(+i k r^2/(2R) - i gouy)
    # (with the e^{ikz} piston dropped, tracked by the path ledger)
    wl, dx = 1e-6, 4e-6
    w, b = _gaussian_wavefront(wl, dx)
    dz = b.rayleigh_distance

    out = lentil.propagate_ptp(w, dz, oversample=2)
    E = out.field
    cc = E.shape[0]//2

    wz, Rz, gouy = b.radius(dz), b.phase_radius(dz), b.gouy(dz)
    k = 2*np.pi/wl

    # on-axis amplitude and Gouy phase
    assert np.isclose(np.abs(E[cc, cc]), b.waist_radius/wz, rtol=1e-9)
    assert np.isclose(np.angle(E[cc, cc]), -gouy, atol=1e-9)

    # radial envelope and quadratic phase curvature
    rr = np.arange(60)*dx
    prof = E[cc, cc:cc+60]
    envelope = (b.waist_radius/wz)*np.exp(-rr**2/wz**2)
    assert np.allclose(np.abs(prof), envelope, rtol=1e-6)
    dphi = np.angle(prof*np.conj(E[cc, cc]))
    assert np.allclose(dphi, np.mod(k*rr**2/(2*Rz) + np.pi, 2*np.pi) - np.pi,
                       atol=1e-6)


@pytest.mark.filterwarnings('ignore:pilot beam diameter')
def test_ptp_talbot():
    # a periodic amplitude grating self-images at the Talbot distance and
    # produces a half-period-shifted image at half the Talbot distance.
    # oversample=1 with a commensurate period makes the FFT wraparound act
    # as an exact periodic boundary (infinite grating)
    wl, dx, n = 1e-6, 10e-6, 256
    period_px = 8
    p = period_px*dx
    zt = 2*p**2/wl

    x = np.arange(n) - n//2
    grating = np.tile(0.5*(1 + np.cos(2*np.pi*x/period_px)), (n, 1))

    w = lentil.Wavefront(wl, pixelscale=dx)
    w = w * lentil.Plane(amplitude=grating, pixelscale=dx)

    i0 = w.intensity
    i_zt = lentil.propagate_ptp(w, zt, oversample=1).intensity
    i_half = lentil.propagate_ptp(w, zt/2, oversample=1).intensity

    assert np.allclose(i_zt, i0, atol=1e-9)
    assert np.allclose(i_half, np.roll(i0, period_px//2, axis=1), atol=1e-9)


def test_ptp_roundtrip():
    w, b = _gaussian_wavefront()
    dz = 0.5*b.rayleigh_distance
    out = lentil.propagate_ptp(lentil.propagate_ptp(w, dz), -dz, oversample=1)
    assert np.allclose(out.field, lentil.pad(w.field, out.shape), atol=1e-10)
    assert out.z == w.z


def test_ptp_energy():
    w, b = _gaussian_wavefront()
    out = lentil.propagate_ptp(w, b.rayleigh_distance, oversample=2)
    assert np.isclose(np.sum(out.intensity), np.sum(w.intensity), rtol=1e-9)


def test_ptp_fft_dft_agree():
    w, b = _gaussian_wavefront(n=128, w0_px=16)
    dz = 0.2*b.rayleigh_distance
    out_fft = lentil.propagate_ptp(w, dz, method='fft', oversample=2)
    out_dft = lentil.propagate_ptp(w, dz, method='dft', oversample=2)
    assert np.allclose(out_fft.field, out_dft.field, atol=1e-10)


@pytest.mark.filterwarnings('ignore:pilot beam diameter')
def test_ptp_state():
    w, b = _gaussian_wavefront(n=128, w0_px=16)
    w.focal_length = 3
    dz = 0.1
    out = lentil.propagate_ptp(w, dz)
    assert out.z == dz
    assert out.path == dz
    assert out.z_focus == 3           # fixed during propagation
    assert out.focal_length == 3 - dz  # derived from current position
    assert out.ptype == lentil.none
    assert np.all(out.pixelscale == w.pixelscale)
    assert out.reference == 'planar'
    assert out.pilot is w.pilot


def test_ptp_accepts_conjugates_returns_none():
    amp = lentil.circle((64, 64), 24)
    w = lentil.Wavefront(1e-6)
    w = w * lentil.Pupil(amplitude=amp, pixelscale=1e-4, focal_length=1)
    assert w.ptype == lentil.pupil
    out = lentil.propagate_ptp(w, 0.01)
    assert out.ptype == lentil.none


def test_ptp_uniform_field_unchanged():
    # a uniform plane wave is invariant under plane-to-plane propagation
    w = lentil.Wavefront(1e-6, pixelscale=1e-5)
    out = lentil.propagate_ptp(w, 0.1)
    assert np.array_equal(out.field, 1+0j)
    assert out.z == 0.1


def test_ptp_invalid():
    w, _ = _gaussian_wavefront(n=64, w0_px=8)
    with pytest.raises(ValueError):
        lentil.propagate_ptp(w, 0.1, method='czt')

    w.reference = 'spherical'
    with pytest.raises(ValueError):
        lentil.propagate_ptp(w, 0.1)

    w2 = lentil.Wavefront(1e-6)  # no pixelscale
    with pytest.raises(ValueError):
        lentil.propagate_ptp(w2, 0.1)


def test_ptp_spread_warning():
    # a fast-diverging pilot beam should trip the guard band diagnostic
    wl, dx, n = 1e-6, 1e-6, 64
    amp = lentil.circle((n, n), 8)
    w = lentil.Wavefront(wl, diameter=16e-6)
    w = w * lentil.Plane(amplitude=amp, pixelscale=dx)
    with pytest.warns(UserWarning, match='pilot beam diameter'):
        lentil.propagate_ptp(w, 0.1)


def test_ptp_tilt_matches_sampled_phasor():
    # end-to-end check of the analytic tilt contract: a plane with tilted
    # OPD propagated with fit_tilt() bookkeeping (integer offset + subpixel
    # kernel phase) matches the same OPD propagated as a sampled phasor.
    # A full-array Gaussian envelope is used so both paths run on
    # identical grids, isolating the tilt handling from wraparound
    wl, dx, n = 1e-6, 4e-6, 256
    dz = 0.05

    f = np.fft.fftshift(np.fft.fftfreq(n, 1/(n*dx)))
    r, c = np.meshgrid(f, f, indexing='ij')
    amp = np.exp(-(r**2 + c**2)/(30*dx)**2)
    # tilt magnitude chosen to give a non-integer ~7.3 px displacement
    theta = 7.3*dx/dz
    opd = theta*c  # c is the x/column coordinate in meters

    p_sampled = lentil.Plane(amplitude=amp, opd=opd, pixelscale=dx)
    p_fit = lentil.Plane(amplitude=amp, opd=opd, pixelscale=dx).fit_tilt()

    w_sampled = lentil.Wavefront(wl) * p_sampled
    w_fit = lentil.Wavefront(wl) * p_fit
    assert w_fit.data[0].tilt  # tilt was actually fit and bookkept

    i_sampled = lentil.propagate_ptp(w_sampled, dz, oversample=2).intensity
    out_fit = lentil.propagate_ptp(w_fit, dz, oversample=2)
    i_fit = np.zeros_like(i_sampled)
    i_fit = out_fit.insert(i_fit)

    assert np.allclose(i_fit, i_sampled, atol=1e-6*np.max(i_sampled))


@pytest.mark.filterwarnings('ignore:pilot beam diameter')
def test_ptp_tilt_persists():
    w, _ = _gaussian_wavefront(n=64, w0_px=8)
    for field in w.data:
        field.tilt.append(lentil.Tilt(x=1e-6, y=0))
    out = lentil.propagate_ptp(w, 0.01)
    assert out.data[0].tilt  # still attached; the beam is still tilted


def test_ptp_per_field_vs_consolidated():
    # segmented, untilted wavefront: per-field (dft) propagation agrees
    # coherently with consolidated (fft) propagation, converging as the
    # per-segment guard band grows (the accuracy/cost trade documented in
    # ADR 0002)
    wl, dx = 1e-6, 4e-6
    segmask = lentil.hex_segments(rings=1, seg_radius=32, seg_gap=2,
                                  flatten=False)
    amp = np.sum(segmask, axis=0)
    p = lentil.Plane(amplitude=amp, mask=segmask, pixelscale=dx)
    w = lentil.Wavefront(wl) * p
    assert len(w.data) == 6

    dz = 1e-3
    out_c = lentil.propagate_ptp(w, dz, method='fft', oversample=8)

    err = {}
    for ovs in (2, 4, 8):
        out_pf = lentil.propagate_ptp(w, dz, method='dft', oversample=ovs)
        i_pf = out_pf.insert(np.zeros(out_c.shape))
        err[ovs] = np.max(np.abs(i_pf - out_c.intensity))/np.max(out_c.intensity)

    assert err[2] < 5e-2
    assert err[4] < err[2]/4
    assert err[8] < err[4]/4


def test_farfield_advances_z():
    amp = lentil.normalize_power(lentil.circle((128, 128), 48))
    w = lentil.Wavefront(1e-6, z=1)
    w = w * lentil.Pupil(amplitude=amp, pixelscale=1e-4, focal_length=2)
    out = lentil.propagate_dft(w, pixelscale=5e-6, shape=32)
    assert out.z == 3
    assert out.path == 2
    assert out.focal_length == 2  # conjugate hop convention
    assert out.reference == 'planar'


def test_defocus_idiom_matches_sampled_defocus():
    # the Q13 seam: far-field to focus then PTP defocus matches an
    # all-far-field model with equivalent sampled defocus OPD
    # (-W020*rho^2 with W020 = dz/(8 F#^2) for PTP by +dz)
    wl, f, D = 1e-6, 1.0, 0.01
    n, rad = 256, 100
    dx = D/(2*rad)
    fno = f/D
    du = wl*fno/4  # Q = 4
    amp = lentil.normalize_power(lentil.circle((n, n), rad))
    dz = 2*wl*fno**2  # W020 = lambda/4

    w = lentil.Wavefront(wl)
    w = w * lentil.Pupil(amplitude=amp, pixelscale=dx, focal_length=f)
    w = lentil.propagate_dft(w, pixelscale=du, shape=64, oversample=2)
    assert w.z == f
    w2 = lentil.propagate_ptp(w, dz, oversample=2)
    assert w2.z == f + dz
    psf_nf = w2.intensity
    c = psf_nf.shape[0]//2
    psf_nf = psf_nf[c-64:c+64, c-64:c+64]

    r, cc = np.meshgrid(np.arange(n)-n//2, np.arange(n)-n//2, indexing='ij')
    rho2 = (r**2 + cc**2)/rad**2
    w020 = dz/(8*fno**2)
    opd = -w020*rho2*(amp > 0)
    wf = lentil.Wavefront(wl)
    wf = wf * lentil.Pupil(amplitude=amp, opd=opd, pixelscale=dx, focal_length=f)
    psf_ff = lentil.propagate_dft(wf, pixelscale=du, shape=64, oversample=2).intensity

    assert np.max(np.abs(psf_nf - psf_ff)) < 2e-3*np.max(psf_ff)


def test_defocus_sign_convention():
    # discriminate the defocus sign with a symmetry-breaking aberration
    # (astigmatism + coma): PTP by +dz matches -W020*rho^2, not +W020*rho^2
    wl, f, D = 1e-6, 1.0, 0.01
    n, rad = 256, 100
    dx = D/(2*rad)
    fno = f/D
    du = wl*fno/4
    amp = lentil.normalize_power(lentil.circle((n, n), rad))
    mask = lentil.circle((n, n), rad) > 0
    ab = 0.2*wl*lentil.zernike(mask, 8) + 0.15*wl*lentil.zernike(mask, 6)
    dz = 2*wl*fno**2

    w = lentil.Wavefront(wl)
    w = w * lentil.Pupil(amplitude=amp, opd=ab, pixelscale=dx, focal_length=f)
    w = lentil.propagate_dft(w, pixelscale=du, shape=64, oversample=2)
    psf_nf = lentil.propagate_ptp(w, dz, oversample=2).intensity
    c = psf_nf.shape[0]//2
    psf_nf = psf_nf[c-64:c+64, c-64:c+64]

    r, cc = np.meshgrid(np.arange(n)-n//2, np.arange(n)-n//2, indexing='ij')
    rho2 = (r**2 + cc**2)/rad**2
    w020 = dz/(8*fno**2)

    err = {}
    for sign in (+1, -1):
        wf = lentil.Wavefront(wl)
        wf = wf * lentil.Pupil(amplitude=amp, opd=ab + sign*w020*rho2*mask,
                               pixelscale=dx, focal_length=f)
        psf = lentil.propagate_dft(wf, pixelscale=du, shape=64, oversample=2).intensity
        err[sign] = np.max(np.abs(psf_nf - psf))/np.max(psf)

    assert err[-1] < 5e-2   # correct sign agrees
    assert err[+1] > 0.3    # wrong sign clearly does not
