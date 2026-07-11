"""Near-field (Fresnel) propagation support.

This module implements the near-field propagation machinery described in
[1], including the pilot (surrogate) Gaussian beam used to select reference
surfaces and propagator composition.

Sign conventions
----------------
Lentil adopts the following conventions for beams traveling in the +z
direction:

* The phase radius of curvature ``R(z) = (z - z_w0) + z_R**2/(z - z_w0)``
  is positive for a diverging beam (waist behind, at smaller z), negative
  for a converging beam (waist ahead, at larger z), and infinite at the
  waist.
* A thin lens with focal length ``f`` (positive = converging) transforms
  the phase radius according to ``1/R2 = 1/R1 - 1/f``.
* A wavefront's ``focal_length`` (distance to the analytic focus) relates
  to the local phase radius by ``focal_length = -R``.

References
----------
[1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
    Engineering v11, ch. 3 (1992)

"""

import warnings

import numpy as np

import lentil
import lentil.field
from lentil.field import Field
import lentil.fourier


def propagate_ptp(wavefront, dz, method='fft', shape=None, oversample=2,
                  scratch=None):
    """Propagate a Wavefront plane-to-plane using the angular spectrum
    method.

    Plane-to-plane (PTP) propagation preserves the wavefront sampling:
    the output pixelscale always equals the input pixelscale. ``oversample``
    pads the transform grid to provide an aliasing guard band around the
    field data [1]; it does not change the output sampling. This differs
    from the far-field propagators, where ``oversample`` refines the output
    sampling.

    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        Wavefront to propagate. Must have a planar reference surface.
    dz : float
        Propagation distance in meters. May be negative (backward
        propagation via the inverse transform).
    method : {'fft', 'dft'}, optional
        Numerical method used to evaluate the transforms. ``'fft'``
        (default) consolidates untilted fields into a single transform;
        ``'dft'`` propagates each field independently on its own
        guard-banded grid, preserving Wavefront sparsity.
    shape : int or (2,) tuple of ints, optional
        Output canvas shape (before oversampling). If None (default), the
        wavefront shape is used. The transform grid is
        ``shape * oversample``.
    oversample : float, optional
        Guard band factor applied to the transform grid. Default is 2,
        which keeps the field at or below 50% of the grid width as
        recommended in [1].
    scratch : complex ndarray, optional
        A pre-allocated array used for zero-padding in the consolidated
        ``'fft'`` path. Providing a sufficiently large scratch array can
        improve broadband propagation performance.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The propagated Wavefront. ptype is always ``lentil.none`` (the
        wavefront is between conjugates); conjugate status is re-asserted
        by multiplying with a Pupil or Image plane on arrival. Tilt
        remains attached to the propagated fields: the geometric
        displacement ``dz*tilt`` is applied via integer Field offsets and
        subpixel linear phase, and the beam remains tilted for subsequent
        propagations.

    References
    ----------
    [1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
        Engineering v11, ch. 3 (1992)

    """
    if method not in ('fft', 'dft'):
        raise ValueError(f"method must be 'fft' or 'dft', got '{method}'")
    if wavefront.reference != 'planar':
        raise ValueError("propagate_ptp requires a planar-reference "
                         "wavefront. Use propagate_fresnel_* to manage "
                         "reference surface transitions.")
    if wavefront.pixelscale is None:
        raise ValueError("wavefront must have a defined pixelscale")
    if wavefront.ptype not in (lentil.none, lentil.pupil, lentil.image):
        raise TypeError("Wavefront must have ptype 'none', 'pupil' "
                        "or 'image'")

    dx = wavefront.pixelscale
    z_out = wavefront.z + dz

    shape = np.asarray(wavefront.shape) if shape is None else np.broadcast_to(shape, (2,))
    shape_out = np.round(np.asarray(shape)*oversample).astype(int)

    # focal_length is derived from z_focus and must be recomputed at the
    # output position; z_focus itself is fixed during propagation
    if wavefront.z_focus is None:
        focal_length_out = None
    else:
        focal_length_out = wavefront.z_focus - z_out

    out = lentil.Wavefront.empty(wavelength=wavefront.wavelength,
                                 pixelscale=dx,
                                 focal_length=focal_length_out,
                                 shape=tuple(shape_out),
                                 ptype=lentil.none,
                                 z=z_out,
                                 pilot=wavefront.pilot,
                                 reference='planar',
                                 path=wavefront.path + dz)

    # sort fields: uniform (scalar) fields pass through PTP unchanged,
    # untilted sampled fields may be consolidated (fft), tilted sampled
    # fields are always propagated per-field to preserve the analytic
    # tilt bookkeeping
    consolidate, per_field = [], []
    for field in wavefront.data:
        if field.size == 1:
            out.data.append(Field(data=field.data, pixelscale=field.pixelscale,
                                  offset=field.offset, tilt=field.tilt))
        elif method == 'fft' and not field.tilt:
            consolidate.append(field)
        else:
            per_field.append(field)

    if consolidate:
        grid_shape = tuple(shape_out)
        if scratch is not None:
            if not all(np.asarray(scratch.shape) >= grid_shape):
                raise ValueError(f'scratch must have shape greater than or '
                                 f'equal to {grid_shape}')
            field_in = scratch[0:grid_shape[0], 0:grid_shape[1]]
            field_in[:] = 0
        else:
            field_in = np.zeros(grid_shape, dtype=complex)
        for field in consolidate:
            field_in = lentil.field.insert(field, field_in)

        data = _ptp(field_in, dx, wavefront.wavelength, dz, method)
        out.data.append(Field(data=data, pixelscale=dx, offset=(0, 0)))

    for field in per_field:
        grid_shape = tuple(np.round(np.asarray(field.shape)*oversample).astype(int))

        # geometric displacement due to tilt: dz*theta, in units of
        # output pixels (PTP output sampling = input sampling, so
        # oversample=1 here)
        shift = field.shift(z=dz, wavelength=wavefront.wavelength,
                            pixelscale=dx, oversample=1, indexing='ij')
        fix_shift = np.fix(shift)
        subpx_shift = shift - fix_shift

        # center the field data in its own guard-banded grid. Note this
        # intentionally uses field.insert (not lentil.pad) so that the
        # placement convention is identical to the consolidated path and
        # to the Field offset bookkeeping for odd-sized data
        field_in = np.zeros(grid_shape, dtype=complex)
        field_in = lentil.field.insert(Field(data=field.data), field_in)
        data = _ptp(field_in, dx, wavefront.wavelength, dz, method,
                    shift=subpx_shift)

        offset = np.asarray(field.offset) + fix_shift.astype(int)
        out.data.append(Field(data=data, pixelscale=dx,
                              offset=tuple(offset), tilt=field.tilt))

    # pilot-based spread diagnostic (Lawrence's guideline: keep the beam
    # at or below ~50% of the array width)
    if out.pilot is not None and shape_out.size > 0:
        width = np.min(shape_out * np.asarray(dx))
        if 2*out.pilot.radius(z_out) > 0.5*width:
            warnings.warn(f'pilot beam diameter '
                          f'{2*out.pilot.radius(z_out):.3g} m exceeds 50% '
                          f'of the propagation grid width {width:.3g} m at '
                          f'z = {z_out:.3g} m; expect aliasing. Increase '
                          f'shape or oversample.')

    return out


def _ptp(field, pixelscale, wavelength, dz, method, shift=(0, 0)):
    # Angular spectrum propagation of a single sampled field:
    # ifft(fft(field) * T(dz)), with any subpixel shift folded into the
    # transfer function as a linear phase
    n = field.shape
    T = _ptp_transfer(n, pixelscale, wavelength, dz, shift)
    if method == 'fft':
        return _ifft2(_fft2(field) * T)
    else:
        alpha = (1/n[0], 1/n[1])
        F = lentil.fourier.dft2(field, alpha=alpha, unitary=True)
        return lentil.fourier.idft2(F * T, alpha=alpha, unitary=True)


def _ptp_transfer(shape, pixelscale, wavelength, dz, shift=(0, 0)):
    # Fresnel angular spectrum transfer function (Lawrence Eq. 87)
    # T(dz) = exp(-j*pi*wavelength*dz*rho**2)
    # with the constant exp(j*k*dz) piston dropped (tracked by the
    # wavefront path ledger). A nonzero shift (in pixels, (r, c)) adds
    # the linear phase that displaces the output field via the Fourier
    # shift theorem.
    pixelscale = np.broadcast_to(pixelscale, (2,))
    fr = np.fft.fftshift(np.fft.fftfreq(shape[0], d=pixelscale[0]))
    fc = np.fft.fftshift(np.fft.fftfreq(shape[1], d=pixelscale[1]))
    FR, FC = np.meshgrid(fr, fc, indexing='ij')
    T = np.exp(-1j*np.pi*wavelength*dz*(FR**2 + FC**2))
    if np.any(np.asarray(shift) != 0):
        dr = shift[0]*pixelscale[0]
        dc = shift[1]*pixelscale[1]
        T = T * np.exp(-2j*np.pi*(FR*dr + FC*dc))
    return T


def _fft2(x):
    return np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(x), norm='ortho'))


def _ifft2(x):
    return np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(x), norm='ortho'))


def propagate_wts(wavefront, dz, method='fft', pixelscale=None, shape=None,
                  oversample=2):
    """Propagate a Wavefront from its waist to a spherical-reference plane
    (Lawrence's WTS operator).

    The wavefront's current position is treated as the waist (the center
    of the output reference sphere). The gentle quadratic phase
    ``q(dz) = exp(j*pi*r^2/(lambda*dz))`` is sampled on the (waist-side)
    input and a single Fourier transform is performed; the output is
    defined relative to a spherical reference surface centered at the
    input plane [1].

    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        Wavefront to propagate. Must have a planar reference surface.
    dz : float
        Propagation distance in meters. May be negative.
    method : {'fft', 'dft'}, optional
        Numerical method. ``'dft'`` supports arbitrary output sampling
        and per-field (sparse) propagation; ``'fft'`` pads to approximate
        a requested ``pixelscale``.
    pixelscale : float or (2,) float, optional
        Physical sampling of the output. If None (default), the natural
        sampling ``lambda*|dz|/(N*dx)`` is used. As with the far-field
        propagators, the actual output sampling is ``pixelscale /
        oversample``.
    shape : int or (2,) tuple of ints, optional
        Output shape (before oversampling). If None (default), the
        wavefront shape is used.
    oversample : float, optional
        Number of times to oversample the output plane. Default is 2.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The propagated Wavefront with a spherical reference surface
        centered at the input plane (``z_focus`` = input position).

    References
    ----------
    [1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
        Engineering v11, ch. 3 (1992)

    """
    return _waist_transform(wavefront, dz, method, pixelscale, shape,
                            oversample, direction='wts')


def propagate_stw(wavefront, dz=None, method='fft', pixelscale=None,
                  shape=None, oversample=2):
    """Propagate a spherical-reference Wavefront to its waist (Lawrence's
    STW operator).

    The propagation must terminate at the center of the wavefront's
    reference sphere (its analytic focus, ``z_focus``). A single Fourier
    transform is performed; the reference sphere exactly absorbs the
    Fresnel kernel's input quadratic phase, and the gentle output
    quadratic phase ``q(dz)`` is sampled into the output array [1]. Note
    this differs from the far-field propagators, which drop the output
    quadratic phase (the Fraunhofer convention): STW retains it so that
    near-field composition is exact.

    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        Wavefront to propagate. Must have a spherical reference surface.
    dz : float, optional
        Propagation distance in meters. If None (default), the distance
        to the wavefront's analytic focus (``z_focus - z``) is used —
        the only legal value; providing any other value raises.
    method, pixelscale, shape, oversample :
        See :func:`propagate_wts`.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The propagated Wavefront at the waist with a planar reference
        surface.

    References
    ----------
    [1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
        Engineering v11, ch. 3 (1992)

    """
    if wavefront.z_focus is None:
        raise ValueError('propagate_stw requires a wavefront with '
                         'spherical reference (defined z_focus)')
    dz_focus = wavefront.z_focus - wavefront.z
    if dz is None:
        dz = dz_focus
    elif not np.isclose(dz, dz_focus):
        raise ValueError(f'STW must terminate at the reference focus: '
                         f'dz = {dz} but z_focus - z = {dz_focus}. Use '
                         f'propagate_fresnel_* for arbitrary distances.')
    return _waist_transform(wavefront, dz, method, pixelscale, shape,
                            oversample, direction='stw')


def _waist_transform(wavefront, dz, method, pixelscale, shape, oversample,
                     direction):
    # Shared implementation of the WTS and STW operators. Both are a
    # single Fourier transform between a planar-reference waist plane and
    # a spherical-reference plane a distance |dz| away; they differ in
    # which side the sampled (always gentle) quadratic phase q(dz) is
    # applied to: the input for WTS, the output for STW.
    if method not in ('fft', 'dft'):
        raise ValueError(f"method must be 'fft' or 'dft', got '{method}'")
    if wavefront.pixelscale is None:
        raise ValueError('wavefront must have a defined pixelscale')
    if dz == 0:
        raise ValueError('dz must be nonzero')
    if direction == 'wts' and wavefront.reference != 'planar':
        raise ValueError('propagate_wts requires a planar-reference '
                         'wavefront')
    if direction == 'stw' and wavefront.reference != 'spherical':
        raise ValueError('propagate_stw requires a spherical-reference '
                         'wavefront')

    dx = wavefront.pixelscale
    wl = wavefront.wavelength
    z_out = wavefront.z + dz
    z_w = wavefront.z if direction == 'wts' else wavefront.z_focus

    shape = np.asarray(wavefront.shape) if shape is None else np.broadcast_to(shape, (2,))
    shape_out = np.round(np.asarray(shape)*oversample).astype(int)

    if pixelscale is None:
        # natural (FFT) output sampling for the oversampled grid
        du = wl*np.abs(dz)/(shape_out*np.asarray(dx))
    else:
        du = np.broadcast_to(pixelscale, (2,))/oversample

    if method == 'fft':
        # the FFT can only produce its natural sampling; pad the
        # transform to (approximately) satisfy the requested output
        # sampling, and recompute the actual sampling from the padded
        # size
        fft_shape = np.round(wl*np.abs(dz)/(np.asarray(dx)*du)).astype(int)
        if np.any(fft_shape < shape_out):
            raise ValueError(f'requested sampling requires transform '
                             f'shape {tuple(fft_shape)} smaller than '
                             f'output shape {tuple(shape_out)}')
        du = wl*np.abs(dz)/(fft_shape*np.asarray(dx))

    # Each waist transform drops the 1/j prefactor of the Fresnel
    # integral (the unitary transform does not include it), leaving the
    # array with an unphysical +pi/2*sign(dz) piston. The path ledger
    # records the correction as its path-equivalent, -sign(dz)*lambda/4,
    # so the smooth Gouy evolution stays in the array while the dropped
    # constant stays recoverable. (Verified against the analytic Gaussian
    # in tests/test_fresnel.py.)
    out = lentil.Wavefront.empty(wavelength=wl,
                                 pixelscale=tuple(du),
                                 focal_length=z_w - z_out,
                                 shape=tuple(shape_out),
                                 ptype=lentil.none,
                                 z=z_out,
                                 pilot=wavefront.pilot,
                                 reference='spherical' if direction == 'wts' else 'planar',
                                 path=wavefront.path + dz - np.sign(dz)*wl/4)
    if direction == 'stw':
        # the output is at the waist; the reference collapses to the
        # current position (focal_length = 0)
        out.focal_length = 0

    for field in wavefront.data:
        if field.size == 1:
            raise ValueError(f'propagate_{direction} requires sampled '
                             f'field data')

        # geometric tilt displacement dz*theta in output pixels
        shift = field.shift(z=dz, wavelength=wl, pixelscale=du,
                            oversample=1, indexing='ij')
        fix_shift = np.fix(shift)
        subpx_shift = shift - fix_shift

        if direction == 'wts':
            data = field.data * _chirp(field.shape, field.offset, dx, wl, dz)
        else:
            data = field.data

        if method == 'dft':
            alpha = (dx[0]*du[0]/(wl*dz), dx[1]*du[1]/(wl*dz))
            F = lentil.fourier.dft2(f=data, alpha=alpha, shape=shape_out,
                                    shift=subpx_shift, offset=field.offset,
                                    unitary=True)
            offset_out = fix_shift.astype(int)
        else:
            grid = np.zeros(tuple(fft_shape), dtype=complex)
            grid = lentil.field.insert(Field(data=data, offset=field.offset), grid)
            if np.any(subpx_shift != 0):
                # displace the output by the subpixel remainder via an
                # input-plane linear phase
                grid *= _linear_phase(tuple(fft_shape), dx, wl, dz,
                                      subpx_shift*du)
            F = _fft2(grid) if dz > 0 else _ifft2(grid)
            offset_out = fix_shift.astype(int)

        if direction == 'stw':
            F = F * _chirp(F.shape, offset_out, du, wl, dz)

        out.data.append(Field(data=F, pixelscale=tuple(du),
                              offset=tuple(offset_out), tilt=field.tilt))

    return out


def _chirp(shape, offset, pixelscale, wavelength, dz):
    # sampled quadratic phase q(dz) = exp(j*pi*r^2/(lambda*dz)) evaluated
    # at grid coordinates offset from the optical axis by `offset` pixels.
    # This is the waist-side chirp of the WTS/STW operators; it is always
    # gentle (~<1 rad across the beam) when |dz| >= z_R
    pixelscale = np.broadcast_to(pixelscale, (2,))
    r = (np.arange(shape[0]) - shape[0]//2 + offset[0])*pixelscale[0]
    c = (np.arange(shape[1]) - shape[1]//2 + offset[1])*pixelscale[1]
    rsq = r[:, np.newaxis]**2 + c[np.newaxis, :]**2
    return np.exp(1j*np.pi*rsq/(wavelength*dz))


def _linear_phase(shape, pixelscale, wavelength, dz, delta):
    # input-plane linear phase that displaces the transform output by
    # delta (meters, (r, c)): exp(+2j*pi*(x1.delta)/(lambda*dz))
    pixelscale = np.broadcast_to(pixelscale, (2,))
    r = (np.arange(shape[0]) - shape[0]//2)*pixelscale[0]
    c = (np.arange(shape[1]) - shape[1]//2)*pixelscale[1]
    ph = r[:, np.newaxis]*delta[0] + c[np.newaxis, :]*delta[1]
    return np.exp(2j*np.pi*ph/(wavelength*dz))


def propagate_fresnel_dft(wavefront, dz, pixelscale=None, shape=None,
                          oversample=2):
    """Propagate a Wavefront a signed distance ``dz`` in the near field,
    automatically composing the PTP/WTS/STW primitive propagators.

    Transforms are evaluated with the matrix DFT, giving independent
    control over output sampling and shape (see
    :func:`~lentil.propagate_dft` for the far-field analog).

    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        Wavefront to propagate. Must carry a pilot beam (see
        :class:`GaussianBeam`); one is created automatically when a
        Wavefront or Plane defines a diameter.
    dz : float
        Signed propagation distance in meters.
    pixelscale : float or (2,) float, optional
        Physical sampling of the output. If None (default), the natural
        sampling of the composed propagation is used. A requested
        pixelscale is honored exactly when the composition contains a
        WTS/STW leg; a pure plane-to-plane propagation cannot change
        sampling and will raise if the request conflicts with the input
        sampling.
    shape : int or (2,) tuple of ints, optional
        Output shape (before oversampling). If None (default), the
        wavefront shape is used.
    oversample : float, optional
        Number of times to oversample the output plane. Default is 2.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`

    See Also
    --------
    propagate_fresnel_fft, propagate_ptp, propagate_wts, propagate_stw
    """
    return _propagate_fresnel(wavefront, dz, 'dft', pixelscale, shape,
                              oversample)


def propagate_fresnel_fft(wavefront, dz, pixelscale=None, shape=None,
                          oversample=2, scratch=None):
    """Propagate a Wavefront a signed distance ``dz`` in the near field,
    automatically composing the PTP/WTS/STW primitive propagators.

    Transforms are evaluated with the FFT. A requested ``pixelscale`` is
    approximately honored by zero-padding (subject to FFT rounding), as
    in :func:`~lentil.propagate_fft`.

    Parameters
    ----------
    wavefront, dz, pixelscale, shape, oversample :
        See :func:`propagate_fresnel_dft`.
    scratch : complex ndarray, optional
        Pre-allocated scratch array passed to the consolidated FFT
        propagation paths.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`

    See Also
    --------
    propagate_fresnel_dft, propagate_ptp, propagate_wts, propagate_stw
    """
    return _propagate_fresnel(wavefront, dz, 'fft', pixelscale, shape,
                              oversample, scratch=scratch)


def _propagate_fresnel(wavefront, dz, method, pixelscale, shape, oversample,
                       scratch=None):
    # Compose the primitive propagators to move a wavefront by dz.
    #
    # Composition is selected by the wavefront's reference surface (which
    # states what the array *means*) and by whether the destination lies
    # inside the pilot beam's Rayleigh distance of the waist:
    #
    #   planar,    destination inside:   PTP(dz)                   [II]
    #   planar,    destination outside:  PTP(z_w-z1), WTS(z2-z_w)  [IO]
    #   spherical, destination inside:   STW(z_w-z1), PTP(z2-z_w)  [OI]
    #   spherical, destination outside:  STW(z_w-z1), WTS(z2-z_w)  [OO]
    #
    # For a spherical-reference wavefront the waist is its analytic focus
    # (z_focus); this differs from Lawrence, who uses the pilot beam's
    # waist. Centering the intermediate plane on the analytic focus makes
    # the reference sphere absorb the wavefront's ideal curvature exactly
    # (nothing is sampled for an ideal beam); the pilot is used only to
    # classify the destination and size guard bands.
    if wavefront.pilot is None:
        raise ValueError('near-field propagation requires a pilot beam. '
                         'Provide a diameter (or explicit pilot) when '
                         'creating the Wavefront, or multiply by a Plane '
                         'that defines one.')
    if wavefront.pixelscale is None:
        raise ValueError('wavefront must have a defined pixelscale')
    if dz == 0:
        raise ValueError('dz must be nonzero')

    z1 = wavefront.z
    z2 = z1 + dz
    zr = wavefront.pilot.rayleigh_distance

    if wavefront.reference == 'spherical':
        if wavefront.z_focus is None:
            raise ValueError("spherical-reference wavefront has no "
                             "z_focus; state is inconsistent")
        z_w = wavefront.z_focus
    else:
        z_w = wavefront.pilot.waist_position

    inside_out = np.abs(z2 - z_w) <= zr

    if wavefront.reference == 'planar':
        if inside_out:
            # II
            if pixelscale is not None and not np.allclose(
                    np.broadcast_to(pixelscale, (2,)), wavefront.pixelscale):
                raise ValueError('plane-to-plane propagation cannot change '
                                 f'sampling: requested pixelscale '
                                 f'{pixelscale} != wavefront pixelscale '
                                 f'{tuple(wavefront.pixelscale)}')
            return propagate_ptp(wavefront, dz, method=method, shape=shape,
                                 oversample=oversample, scratch=scratch)
        else:
            # IO
            if not np.isclose(z_w, z1):
                wavefront = propagate_ptp(wavefront, z_w - z1, method=method,
                                          oversample=oversample,
                                          scratch=scratch)
            return propagate_wts(wavefront, z2 - z_w, method=method,
                                 pixelscale=pixelscale, shape=shape,
                                 oversample=oversample)
    else:
        if inside_out:
            # OI: the STW leg lands the requested sampling, which the
            # final PTP leg preserves
            wavefront = propagate_stw(wavefront, method=method,
                                      pixelscale=pixelscale, shape=shape,
                                      oversample=oversample)
            if np.isclose(z2, z_w):
                return wavefront
            # guard band was established by the STW leg
            return propagate_ptp(wavefront, z2 - z_w, method=method,
                                 shape=wavefront.shape, oversample=1,
                                 scratch=scratch)
        else:
            # OO
            wavefront = propagate_stw(wavefront, method=method,
                                      oversample=oversample)
            return propagate_wts(wavefront, z2 - z_w, method=method,
                                 pixelscale=pixelscale, shape=shape,
                                 oversample=oversample)


class GaussianBeam:
    """Gaussian beam parameters.

    Primarily used as the *pilot beam* carried by a
    :class:`~lentil.Wavefront` (as ``wavefront.pilot``) to steer near-field
    propagation decisions: reference surface selection and propagator
    composition. In that role the pilot is a bookkeeping surrogate for the
    actual beam, not physics — it only needs to be representative enough to
    pick sampling regimes sensibly [1].

    Parameters
    ----------
    waist_radius : float
        Beam waist radius (1/e field amplitude radius) in meters
    wavelength : float
        Wavelength in meters
    waist_position : float, optional
        Axial location of the beam waist in meters. Default is 0.

    References
    ----------
    [1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
        Engineering v11, ch. 3 (1992)

    """
    __slots__ = ('waist_radius', 'wavelength', 'waist_position')

    def __init__(self, waist_radius, wavelength, waist_position=0):
        if waist_radius <= 0:
            raise ValueError('waist_radius must be positive')
        if wavelength <= 0:
            raise ValueError('wavelength must be positive')

        #: float: Beam waist radius
        self.waist_radius = waist_radius

        #: float: Wavelength
        self.wavelength = wavelength

        #: float: Axial location of the beam waist
        self.waist_position = waist_position

    def __repr__(self):
        return (f'{self.__class__.__name__}(waist_radius='
                f'{self.waist_radius:.4g}, wavelength='
                f'{self.wavelength:.4g}, waist_position='
                f'{self.waist_position:.4g})')

    @property
    def rayleigh_distance(self):
        """Rayleigh distance :math:`z_R = \\pi w_0^2/\\lambda`

        Returns
        -------
        float
        """
        return np.pi * self.waist_radius**2 / self.wavelength

    @property
    def divergence(self):
        """Far-field divergence half-angle :math:`\\theta = w_0/z_R`
        in radians

        Returns
        -------
        float
        """
        return self.waist_radius / self.rayleigh_distance

    def radius(self, z):
        """Beam radius :math:`w(z) = w_0\\sqrt{1 + ((z-z_{w0})/z_R)^2}`

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        float
        """
        return self.waist_radius * np.sqrt(1 + ((z - self.waist_position)/self.rayleigh_distance)**2)

    def phase_radius(self, z):
        """Phase radius of curvature
        :math:`R(z) = (z-z_{w0}) + z_R^2/(z-z_{w0})`

        ``R > 0`` for a diverging beam (waist behind), ``R < 0`` for a
        converging beam (waist ahead). Returns ``inf`` at the waist.

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        float
        """
        dz = z - self.waist_position
        if dz == 0:
            return np.inf
        return dz + self.rayleigh_distance**2/dz

    def gouy(self, z):
        """Gouy phase shift :math:`\\theta(z) = \\arctan((z-z_{w0})/z_R)`
        in radians

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        float
        """
        return np.arctan((z - self.waist_position)/self.rayleigh_distance)

    def inside(self, z):
        """True if ``z`` is inside the Rayleigh distance of the waist

        A position inside the Rayleigh distance uses a planar reference
        surface; a position outside uses a spherical reference surface.

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        bool
        """
        return np.abs(z - self.waist_position) <= self.rayleigh_distance

    @classmethod
    def from_radius(cls, radius, phase_radius, wavelength, z=0):
        """Create a GaussianBeam from beam radius and phase radius known
        at some axial position (waist finding).

        Parameters
        ----------
        radius : float
            Beam radius :math:`w` at ``z``
        phase_radius : float
            Phase radius of curvature :math:`R` at ``z``. May be ``inf``
            (the beam has its waist at ``z``).
        wavelength : float
            Wavelength
        z : float, optional
            Axial position where ``radius`` and ``phase_radius`` are
            defined. Default is 0.

        Returns
        -------
        :class:`GaussianBeam`
        """
        if np.isinf(phase_radius):
            return cls(radius, wavelength, waist_position=z)

        # Lawrence Eq. 56, rewritten in Lentil's sign convention:
        #   z - z_w0 = R/(1 + (lambda*R/(pi*w**2))**2)
        #   w_0 = w/sqrt(1 + (pi*w**2/(lambda*R))**2)
        a = wavelength*phase_radius/(np.pi*radius**2)
        dz = phase_radius/(1 + a**2)
        waist_radius = radius/np.sqrt(1 + 1/a**2)
        return cls(waist_radius, wavelength, waist_position=z-dz)

    def lens(self, focal_length, z):
        """Return a new GaussianBeam transformed by a thin lens at ``z``.

        The beam radius is unchanged by the lens; the phase radius is
        transformed according to ``1/R2 = 1/R1 - 1/f`` (Lawrence Eq. 57 in
        Lentil's sign convention) and the new waist is found from the
        transformed parameters.

        Parameters
        ----------
        focal_length : float
            Lens focal length. Positive is converging. May be ``inf``
            (no optical power).
        z : float
            Axial position of the lens

        Returns
        -------
        :class:`GaussianBeam`
        """
        if np.isinf(focal_length):
            return GaussianBeam(self.waist_radius, self.wavelength,
                                waist_position=self.waist_position)

        w = self.radius(z)
        r1 = self.phase_radius(z)

        if np.isinf(r1):
            r2 = -focal_length
        elif r1 == focal_length:
            # collimating lens - output waist is at the lens
            return GaussianBeam(w, self.wavelength, waist_position=z)
        else:
            r2 = r1*focal_length/(focal_length - r1)

        return GaussianBeam.from_radius(w, r2, self.wavelength, z=z)
