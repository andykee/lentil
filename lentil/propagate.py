import numpy as np

import lentil
from lentil.field import Field
from lentil.wavefront import Wavefront

def propagate_dft(wavefront, shape, pixelscale, prop_shape=None, 
                  oversample=2, inplace=True):
    """Propagate a Wavefront using Fraunhofer diffraction.

    Parameters
    ----------
    shape : int or (2,) tuple of ints
        Shape of output Wavefront.
    pixelscale : float or (2,) float
        Physical sampling of output Wavefront. If a single value is supplied,
        the output is assumed to be uniformly sampled in both x and y.
    prop_shape : int or (2,) tuple of ints, optional
        Shape of propagation output. If None (default),
        ``prop_shape = prop``. If ``prop_shape != prop``, the propagation
        result is placed in the appropriate location in the output plane.
        ``prop_shape`` should not be larger than ``prop``.
    oversample : int, optional
        Number of times to oversample the output plane. Default is 2.
    inplace : bool, optional
        If True (default) the Wavefront is propagated in-place, otherwise
        a copy is created and propagated.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        A Wavefront propagated to the specified image plane
    """
    
    ptype_out = _propagate_ptype(wavefront.ptype, method='fraunhofer')
    
    shape = np.broadcast_to(shape, (2,))
    prop_shape = shape if prop_shape is None else np.broadcast_to(prop_shape, (2,))
    shape_out = shape * oversample
    prop_shape_out = prop_shape * oversample

    dx = wavefront.pixelscale
    du = np.broadcast_to(pixelscale, (2,))
    z = wavefront.focal_length

    data = wavefront.data

    if inplace:
        out = wavefront
        out.data = []
        out.pixelscale = du/oversample
        out.shape = shape_out
        out.ptype = ptype_out
    else:
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
                                       npix=prop_shape_out,
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