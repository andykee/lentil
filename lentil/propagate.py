import numpy as np

import lentil
import lentil.extent
from lentil.field import Field
from lentil.wavefront import Wavefront


def propagate_fft(wavefront, pixelscale, shape=None, oversample=2, 
                  scratch=None):
    """Propagate a Wavefront in the far-field using the FFT.
    
    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        Wavefront to propagate
    pixelscale : float or (2,) float
        Physical sampling of output Wavefront. If a single value is supplied,
        the output is assumed to be uniformly sampled in both x and y.
    shape : int or (2,) tuple of ints or None
        Shape of output Wavefront. If None (default), the wavefront shape is
        used.
    oversample : float, optional
        Number of times to oversample the output plane. Default is 2.
    scratch : complex ndarray, optional
        A pre-allocated array used for padding. Providing a sufficiently
        large scratch array can improve broadband propagation performance by 
        avoiding the repeated allocation of arrays of zeros.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The propagated Wavefront

    """
    if _has_tilt(wavefront):
        raise NotImplementedError('propagate_fft does not support Wavefronts '
                                  'with fitted tilt. Use propagate_dft instead.')

    ptype_out = _propagate_ptype(wavefront.ptype, method='fraunhofer')
    pixelscale = np.broadcast_to(pixelscale, (2,))
    fft_shape, prop_wavelength = _fft_shape(wavefront.pixelscale, 
                                            pixelscale, 
                                            wavefront.focal_length, 
                                            wavefront.wavelength, 
                                            oversample)
    
    # TODO: verify field_shape is smaller than fft_shape
    # TODO: what if field is bigger than fft_shape?
    # field_shape = wavefront.shape

    if shape is None:
        shape_out = tuple(fft_shape)
        shape = (fft_shape[0]//oversample, fft_shape[1]//oversample)
    else:
        shape = tuple(np.broadcast_to(shape, (2,)))
        if np.any(shape > fft_shape/oversample):
            raise ValueError(f'requested shape {tuple(shape)} is larger in at '
                            f'least one dimension than maximum propagation '
                            f'shape {tuple(fft_shape//oversample)}')
        else:
            shape_out = (shape[0] * oversample, shape[1]*oversample)

    out = Wavefront.empty(wavelength=prop_wavelength,
                          pixelscale = pixelscale/oversample,
                          focal_length=wavefront.focal_length,
                          shape = shape_out,
                          ptype = ptype_out)
    
    if scratch is not None:
        if not all(np.asarray(scratch.shape) > fft_shape):
             raise ValueError(f'scratch must have shape greater than or '
                              f'equal to {tuple(fft_shape)}')

        # zero out the portion of scratch that we're going to use for the
        # propagation and then insert the Wavefront field(s) into scratch
        scratch[0:fft_shape[0], 0:fft_shape[1]] = 0
        for field in wavefront.data:
            scratch[0:fft_shape[0], 0:fft_shape[1]] = lentil.field.insert(field, scratch[0:fft_shape[0], 0:fft_shape[1]])
        field =_fft2(scratch[0:fft_shape[0], 0:fft_shape[1]])

    else:
        field = lentil.pad(wavefront.field, fft_shape)
        field = _fft2(field)
    
    out.data.append(Field(data=field, pixelscale=pixelscale/oversample))

    return out


def scratch_shape(wavelength, dx, du, z, oversample):
    """Compute the scratch shape required for an FFT propagation
    
    Parameters
    ----------
    wavelength : float or array_like
        Wavelength or list of wavelengths
    dx : float or (2,) float
        Physical sampling of input Plane. If a single value is supplied,
        the input is assumed to be uniformly sampled in both x and y.
    du : float or (2,) float
        Physical sampling of output Wavefront. If a single value is supplied,
        the output is assumed to be uniformly sampled in both x and y.
    z : float
        Propagation distance
    oversample : float
        Number of times to oversample the output plane.

    Returns
    -------
    shape : tuple
        Scratch shape

    """
    dx = np.broadcast_to(dx, (2,))
    du = np.broadcast_to(du, (2,))
    fft_shape, _ = _fft_shape(dx, du, z, np.max(wavelength), oversample)
    return tuple(fft_shape)


def _dft_alpha(dx, du, wavelength, z, oversample):
    return ((dx[0]*du[0])/(wavelength*z*oversample),
            (dx[1]*du[1])/(wavelength*z*oversample))


def _fft_shape(dx, du, z, wavelength, oversample):
    # Compute pad shape to satisfy requested sampling. Propagation wavelength
    # is recomputed to account for integer padding of input plane
    alpha = _dft_alpha(dx, du, z, wavelength, oversample)
    fft_shape = np.round(np.reciprocal(alpha)).astype(int)       
    prop_wavelength = np.min((fft_shape/oversample * dx * du)/z)
    return fft_shape, prop_wavelength


def _fft2(x):
    return np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(x), norm='ortho'))


def _has_tilt(wavefront):
    # Return True if and Wavefront Field has nonempty tilt
    for field in wavefront.data:
        if field.tilt:
            return True
    return False


def propagate_dft(wavefront, pixelscale, shape=None, prop_shape=None, 
                  oversample=2, mask=None):
    """Propagate a Wavefront in the far-field using the DFT.

    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        Wavefront to propagate
    pixelscale : float or (2,) float
        Physical sampling of output Wavefront. If a single value is supplied,
        the output is assumed to be uniformly sampled in both x and y.
    shape : int or (2,) tuple of ints or None
        Shape of output Wavefront. If None (default), the wavefront shape is
        used.
    prop_shape : int or (2,) tuple of ints, optional
        Shape of propagation output. If None (default),
        ``prop_shape = prop``. If ``prop_shape != prop``, the propagation
        result is placed in the appropriate location in the output plane.
        ``prop_shape`` should not be larger than ``prop``.
    oversample : float, optional
        Number of times to oversample the output plane. Default is 2.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The propagated Wavefront
    """
    
    ptype_out = _propagate_ptype(wavefront.ptype, method='fraunhofer')
    
    shape = np.asarray(wavefront.shape) if shape is None else np.broadcast_to(shape, (2,))
    prop_shape = np.asarray(shape) if prop_shape is None else np.broadcast_to(prop_shape, (2,))
    shape_out = shape * oversample
    prop_shape_out = prop_shape * oversample

    if mask is not None:
        mask = np.asarray(mask)
        if np.all(mask.shape != shape_out):
            raise ValueError(f'shape mismatch: mask shape {mask.shape} != output shape {tuple(shape_out)}')
        mask_shape = _mask_shape(mask, threshold=0)
        mask_shift = _mask_shift(mask, threshold=0)
        out_extent = lentil.extent.array_extent(mask_shape, mask_shift)
    else:
        out_extent = lentil.extent.array_extent(shape_out, shift=(0,0))

    dx = wavefront.pixelscale
    du = np.broadcast_to(pixelscale, (2,))
    z = wavefront.focal_length

    data = wavefront.data

    out = Wavefront.empty(wavelength=wavefront.wavelength,
                          pixelscale = du/oversample,
                          focal_length=wavefront.focal_length,
                          shape = shape_out,
                          ptype = ptype_out)
        
    for field in data:
        # compute the field shift from any embedded tilts. note the return value
        # is specified as (r, c)
        shift = field.shift(z=wavefront.focal_length, wavelength=wavefront.wavelength,
                            pixelscale=du, oversample=oversample,
                            indexing='ij')
        
        fix_shift = np.fix(shift)
        subpx_shift = shift - fix_shift

        prop_extent = lentil.extent.array_extent(prop_shape_out, fix_shift)

        if lentil.extent.intersect(out_extent, prop_extent):
            intersect_shape = lentil.extent.intersection_shape(out_extent, prop_extent)
            intersect_shift = lentil.extent.intersection_shift(out_extent, prop_extent)
            intersect_extent = lentil.extent.array_extent(intersect_shape, intersect_shift)

            # compute additional shift rquired to offset any output clipping that
            # may have occurred due to a mask or extending beyond the full output
            # shape. 
            # NOTE: It appears the same functionality is available using 
            # extent.intersection_shift() but there is an off by one error on only
            # one of the shift dimensions for some reason that causes the tests to 
            # fail
            prop_center = lentil.extent.array_center(prop_extent)
            intersect_center = lentil.extent.array_center(intersect_extent)
            prop_shift = np.array(prop_center) - np.array(intersect_center)

            alpha = _dft_alpha(dx=dx, du=du, z=z,
                               wavelength=wavefront.wavelength,
                               oversample=oversample)
            data = lentil.fourier.dft2(f=field.data, alpha=alpha,
                                       shape=intersect_shape,
                                       shift=prop_shift + subpx_shift,
                                       offset=field.offset, unitary=True)
            out.data.append(Field(data=data, pixelscale=du / oversample,
                                  offset=intersect_shift))

    return out


def _mask_shape(x, threshold=0):
    # compute the shape of a masked area inside an array of zeros
    rmin, rmax, cmin, cmax = lentil.boundary(x, threshold)
    return rmax - rmin + 1, cmax - cmin + 1


def _mask_shift(x, threshold=0):
    # compute the shift of a masked area inside an array of zeros
    shape_full = x.shape
    rc_full, cc_full = shape_full[0]//2, shape_full[1]//2

    rmin_extent, rmax_extent, cmin_extent, cmax_extent = lentil.boundary(x, threshold)
    shape_extent = (rmax_extent-rmin_extent+1, cmax_extent-cmin_extent+1)
    rc_extent, cc_extent = rmin_extent + shape_extent[0]//2, cmin_extent + shape_extent[1]//2

    return rc_extent - rc_full, cc_extent - cc_full


def _propagate_ptype(ptype, method='fraunhofer'):
    if method == 'fraunhofer':
        if ptype not in (lentil.pupil, lentil.image):
            raise TypeError("Wavefront must have ptype 'pupil' "\
                        "or 'image'")
        
        if ptype == lentil.pupil:
            return lentil.image
        else:
            return lentil.pupil