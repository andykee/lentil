import numpy as np

import lentil
from lentil.field import Field
from lentil.wavefront import Wavefront

def propagate_dft(wavefront, pixelscale, shape=None, prop_shape=None, 
                  oversample=2, inplace=True):
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
    inplace : bool, optional
        If True (default) the Wavefront is propagated in-place, otherwise
        a copy is created and propagated.

    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        The orioagated Wavefront
    """
    
    ptype_out = _propagate_ptype(wavefront.ptype, method='fraunhofer')
    
    shape = wavefront.shape if shape is None else np.broadcast_to(shape, (2,))
    prop_shape = shape if prop_shape is None else np.broadcast_to(prop_shape, (2,))
    shape_out = np.asarray(shape) * oversample
    prop_shape_out = np.asarray(prop_shape) * oversample

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

        #if _overlap(prop_shape_out, fix_shift, shape_out):
        prop_shape_out, fix_shift = _update_shape_offset(out_shape=shape_out, 
                                                         field_shape=prop_shape_out, 
                                                         field_offset=fix_shift)
        if any(prop_shape_out):
            alpha = lentil.helper.dft_alpha(dx=dx, du=du, z=z,
                                            wave=wavefront.wavelength,
                                            oversample=oversample)
            data = lentil.fourier.dft2(f=field.data, alpha=alpha,
                                       shape=prop_shape_out,
                                       shift=subpx_shift,
                                       cin=field.offset,
                                       cout=(0,0),
                                       unitary=True)
            
            out.data.append(Field(data=data, pixelscale=du/oversample,
                                  offset=fix_shift))
    
    if not out.data:
        out.data.append(Field(data=0))

    return out


def _propagate_ptype(ptype, method='fraunhofer'):
    if method == 'fraunhofer':
        if ptype not in (lentil.pupil, lentil.image):
            raise TypeError("Wavefront must have ptype 'pupil' "\
                        "or 'image'")
        
        if ptype == lentil.pupil:
            return lentil.image
        else:
            return lentil.pupil
        

def _update_shape_offset(out_shape, field_shape, field_offset, mask=None):
    out_shape = tuple(int(n) for n in out_shape)
    field_shape = tuple(int(n) for n in field_shape)
    field_offset = tuple(int(n) for n in field_offset)
    field_ul = (out_shape[0]//2 - field_shape[0]//2 + field_offset[0],
                out_shape[1]//2 - field_shape[1]//2 + field_offset[1])
    
    if _overlap(out_shape, field_shape, field_ul):

        # Field slice indices
        field_rmin = int(0)
        field_rmax = int(field_shape[0])
        field_cmin = int(0)
        field_cmax = int(field_shape[1])

        # Output insertion slice indices
        out_rmin = int(field_ul[0])
        out_rmax = int(field_ul[0] + field_shape[0])
        out_cmin = int(field_ul[1])
        out_cmax = int(field_ul[1] + field_shape[1])

        # reconcile the field and output insertion indices
        if out_rmin < 0:
            field_rmin = -1 * out_rmin
            out_rmin = 0

        if out_rmax > out_shape[0]:
            field_rmax -= out_rmax - out_shape[0]
            out_rmax = out_shape[0]

        if out_cmin < 0:
            field_cmin = -1 * out_cmin
            out_cmin = 0

        if out_cmax > out_shape[1]:
            field_cmax -= out_cmax - out_shape[1]
            out_cmax = out_shape[1]

        out_center = (out_shape[0]//2, out_shape[1]//2)
        
        field_shape = (field_rmax-field_rmin, field_cmax-field_cmin)
        #print('field shape:', field_shape)

        field_center = (out_rmin + field_shape[0]//2,
                        out_cmin + field_shape[1]//2)
        #print('field center:', field_center)

        field_offset = (field_center[0] - out_center[0],
                        field_center[1] - out_center[1])
        #print('field offset:', field_offset)

        return field_shape, field_offset
        
    else:
        return (), field_offset


def _overlap(out_shape, field_shape, field_ul):
    # Return True if there's any overlap between a shifted field and the
    # output shape
    if field_ul[0] > out_shape[0]:
        return False
    if field_ul[0] + field_shape[0] < 0:
        return False
    if field_ul[1] > out_shape[1]:
        return False
    if field_ul[1] + field_shape[1] < 0:
        return False
    return True