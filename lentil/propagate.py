import numpy as np

import lentil
from lentil.field import Field
from lentil.wavefront import Wavefront


def propagate_fft(wavefront, pixelscale, shape=None, oversample=2, 
                  scratch=None):

    alpha = (wavefront.pixelscale * pixelscale)/(wavefront.wavelength * wavefront.focal_length)
    pad_shape = np.round(oversample * np.reciprocal(alpha)).astype(int)

    # since we can only pad to integer precision, the actual wavelength 
    # represented by the propagation may be slightly different
    wavelength_prop = (pad_shape/oversample * wavefront.pixelscale * pixelscale)/wavefront.focal_length

    if shape is None:
        shape_out = pad_shape
    else:
        shape = np.broadcast_to(shape, (2,))
        shape_out = (np.array(shape) * oversample).astype(int)

    ptype_out = _propagate_ptype(wavefront.ptype, method='fraunhofer')

    out = Wavefront.empty(wavelength=wavelength_prop,
                          pixelscale = pixelscale/oversample,
                          shape = shape_out,
                          ptype = ptype_out)

    # if scratch is not None:
    #     if scratch.dtype != complex:
    #         raise TypeError('scratch must be complex')
        
    #     if not all(np.asarray(scratch.shape) > pad_shape):
    #         raise ValueError(f'scratch must have shape >= {pad_shape}')

    #     field = scratch  # field is just a reference to scratch
    #     field[:] = 0

    
    # else:
    #     field = wavefront.field

    # if wavefront.tilt:
    #     raise ValueError

    field = lentil.pad(wavefront.field, pad_shape)
    field = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(field), norm='ortho'))

    if not (field.shape == shape_out).all():
        field = lentil.pad(field, shape_out)
    
    out.data.append(Field(data=field, pixelscale=pixelscale/oversample))

    return out



    


def propagate_dft(wavefront, pixelscale, shape=None, prop_shape=None, 
                  oversample=2):
    """Propagate a Wavefront using Fraunhofer diffraction.

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
    oversample : int, optional
        Number of times to oversample the output plane. Default is 2.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The orioagated Wavefront
    """
    
    ptype_out = _propagate_ptype(wavefront.ptype, method='fraunhofer')
    
    shape = wavefront.shape if shape is None else np.broadcast_to(shape, (2,))
    prop_shape = shape if prop_shape is None else np.broadcast_to(prop_shape, (2,))
    shape_out = shape * oversample
    prop_shape_out = prop_shape * oversample

    dx = wavefront.pixelscale
    du = np.broadcast_to(pixelscale, (2,))
    z = wavefront.focal_length

    data = wavefront.data

    out = Wavefront.empty(wavelength=wavefront.wavelength,
                          pixelscale = du/oversample,
                          shape = shape_out,
                          ptype = ptype_out)
        
    for field in data:
        # compute the field shift from any embedded tilts. note the return value
        # is specified in terms of (r, c)
        shift = field.shift(z=wavefront.focal_length, wavelength=wavefront.wavelength,
                            pixelscale=du, oversample=oversample,
                            indexing='ij')
        
        fix_shift = np.fix(shift)
        subpx_shift = shift - fix_shift

        if _overlap(prop_shape_out, fix_shift, shape_out):
            alpha = lentil.helper.dft_alpha(dx=dx, du=du, z=z,
                                            wave=wavefront.wavelength,
                                            oversample=oversample)
            data = lentil.fourier.dft2(f=field.data, alpha=alpha,
                                       shape=prop_shape_out,
                                       shift=subpx_shift,
                                       offset=field.offset, unitary=True)
            out.data.append(Field(data=data, pixelscale=du/oversample,
                                  offset=fix_shift))
    
    if not out.data:
        out.data.append(Field(data=0))

    return out

def _overlap(field_shape, field_shift, output_shape):
    # Return True if there's any overlap between a shifted field and the
    # output shape
    output_shape = np.asarray(output_shape)
    field_shape = np.asarray(field_shape)
    field_shift = np.asarray(field_shift)

    # Output coordinates of the upper left corner of the shifted data array
    field_shifted_ul = (output_shape / 2) - (field_shape / 2) + field_shift

    if field_shifted_ul[0] > output_shape[0]:
        return False
    if field_shifted_ul[0] + field_shape[0] < 0:
        return False
    if field_shifted_ul[1] > output_shape[1]:
        return False
    if field_shifted_ul[1] + field_shape[1] < 0:
        return False
    return True

def _propagate_ptype(ptype, method='fraunhofer'):
    if method == 'fraunhofer':
        if ptype not in (lentil.pupil, lentil.image):
            raise TypeError("Wavefront must have ptype 'pupil' "\
                        "or 'image'")
        
        if ptype == lentil.pupil:
            return lentil.image
        else:
            return lentil.pupil