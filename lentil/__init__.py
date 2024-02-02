# Lentil
# Heart-healthy physical optics
#
# Copyright (c) 2020-2024, California Institute of Technology ("Caltech"). 
# U.S. Government sponsorship acknowledged.

__version__ = '0.8.0'

from lentil.ptype import ptype

none = ptype('none')
pupil = ptype('pupil')
image = ptype('image')
tilt = ptype('tilt')
transform = ptype('transform')

from lentil.convolvable import jitter, smear

from lentil import detector

from lentil import fourier

from lentil.plane import (
	Plane,
	Pupil,
	Image,
    Tilt,
	DispersiveTilt,
	Grism,
	LensletArray,
	Rotate,
	Flip
)

from lentil.propagate import (
    propagate_dft, 
    propagate_fft,
    scratch_shape
)

from lentil import radiometry

from lentil.segmented import hex_segments

from lentil.shape import (
    circle,
	hexagon,
	rectangle,
    spider
)

from lentil import util
from lentil.util import (
	centroid,
	pad,
	window,
	boundary,
	rebin,
	rescale,
	pixelscale_nyquist,
	min_sampling,
	normalize_power,
	sanitize_shape,
	sanitize_bandpass
)

from lentil.wavefront import Wavefront

from lentil.wfe import power_spectrum, translation_defocus

from lentil.zernike import (
	zernike,
	zernike_compose,
	zernike_fit,
	zernike_remove,
	zernike_basis,
	zernike_coordinates
)
