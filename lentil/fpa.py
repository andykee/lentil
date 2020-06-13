import numpy as np
from scipy import ndimage

from lentil.convolvable import Pixel
from lentil.plane import Detector
from lentil.radiometry import Spectrum
from lentil import detector
from lentil import util

__all__ = ['FPA', 'BayerFPA']


class FPA(Detector):
    """Focal plane array base class.

    The entire signal flow through the FPA model is depicted in the figure
    below. Boxes represented by attributes which are ``None`` are simply passed
    through without modifying the signal in any way.

    .. image:: ../_static/img/detector_model.png

    Parameters
    ----------
    pixelscale : float
        Pixel size in meters. Pixels are assumed to be square.

    shape : int or array_like
        Number of pixels as (rows, cols). If a single value is provided, the
        :class:`FPA` is assumed to be square with nrows = ncols = shape.

    **kwargs : :class:`FPA` properties, optional
        Additional detector attributes which can be optionally specified. Any
        attributes not defined will default to None and have no effect on
        image simulations. The full list of supported attributes is:

        ======================== ======================================================================
        :attr:`charge_diffusion` :class:`~detector.ChargeDiffusion` model
        :attr:`cosmic_rays`      :class:`~detector.CosmicRay` model
        :attr:`dark_shot_noise`  Dark current :class:`~detector.ShotNoise` model
        :attr:`dark_signal`      :class:`~detector.DarkSignal` model
        :attr:`defects`          :class:`~detector.Defect` or list of :class:`~detector.Defect` objects
        :attr:`dsnu`             Dark signal nonuniformity :class:`~detector.FPN` model
        :attr:`gain`             :class:`~detector.Gain` model
        :attr:`offset_fpn`       Offset :class:`~detector.FPN` model
        :attr:`prnu`             Light signal nonuniformity :class:`~detector.FPN` model
        :attr:`qe`               QE as a :class:`~lentil.radiometry.Spectrum` or a scalar
        :attr:`read_noise`       :class:`~detector.ReadNoise` model
        :attr:`shot_noise`       :class:`~detector.ShotNoise` model
        ======================== ======================================================================

    Examples
    --------
    Additional arguments can be passed in individually:

    .. code:: pycon

        >>> f = FPA(pixelscale=5e-6, shape=(256, 256), read_noise=ReadNoise(50))

    They can also be packed into a dictionary and passed in that way:

    .. code:: pycon

        >>> f_attrs = {'shot_noise': ShotNoise(), 'dark_current': DarkCurrent(500)}
        >>> f = FPA(pixelscale=5e-6, shape=(256, 256), **f_attrs)

    """
    # TODO: rewrite this like Plane constructor. Also need to define __init_subclass__
    def __init__(self, pixelscale, shape, **kwargs):
        super().__init__(pixelscale, shape)

        self.qe = kwargs.get('qe', None)
        self.shot_noise = kwargs.get('shot_noise', None)
        self.prnu = kwargs.get('prnu', None)
        self.dark_signal = kwargs.get('dark_signal', None)
        self.dark_shot_noise = kwargs.get('dark_shot_noise', None)
        self.dsnu = kwargs.get('dsnu', None)
        self._gain = kwargs.get('gain', None)
        self._read_noise = kwargs.get('read_noise', None)
        self._offset_fpn = kwargs.get('offset_fpn', None)
        self._charge_diffusion = kwargs.get('charge_diffusion', None)
        self._defects = kwargs.get('defects', None)
        self._cosmic_rays = kwargs.get('cosmic_rays', None)

    @property
    def qe(self):
        """Quantum efficiency curve used to convert detected photons to
        electrons.

        Returns
        -------
        qe :  :class:`~lentil.radiometry.Spectrum`

        """
        return self._qe

    @qe.setter
    def qe(self, value):
        if value is None:
            self._qe = None
        elif isinstance(value, Spectrum):
            self._qe = value
        else:
            raise AttributeError('qe must be a Spectrum object or None')

    @property
    def shot_noise(self):
        """Shot noise model applied to the light signal. If ``None`` (default),
        shot noise is not applied.

        Returns
        -------
        shot_noise : :class:`~lentil.detector.ShotNoise` or None, optional

        """
        return self._shot_noise

    @shot_noise.setter
    def shot_noise(self, value):
        if value is None:
            self._shot_noise = None
        elif isinstance(value, detector.ShotNoise):
            self._shot_noise = value
        else:
            raise AttributeError('shot_noise must be a ShotNoise object or None')

    @property
    def prnu(self):
        """Photo response non-uniformity model. If ``None`` (default), PRNU is
        not applied.

        Returns
        -------
        prnu : :class:`~lentil.detector.FPN` or None, optional

        """
        return self._prnu

    @prnu.setter
    def prnu(self, value):
        if value is None:
            self._prnu = None
        elif isinstance(value, detector.FPN):
            self._prnu = value
        else:
            raise AttributeError('prnu must be an FPN object or None')

    @property
    def dark_signal(self):
        """Dark signal model. If not ``None`` (default), dark current is not
        applied.

        Returns
        -------
        dark_signal : :class:`~lentil.detector.DarkSignal` or None, optional

        """
        return self._dark_signal

    @dark_signal.setter
    def dark_signal(self, value):
        if value is None:
            self._dark_signal = None
        elif isinstance(value, detector.DarkSignal):
            self._dark_signal = value
        else:
            raise AttributeError('dark_signal must be a DarkSignal object or None')

    @property
    def dark_shot_noise(self):
        """Shot noise model applied to the dark signal. If ``None`` (default),
        shot noise is not applied.

        Returns
        -------
        dark_shot_noise : :class:`~lentil.detector.ShotNoise` or None, optional

        """
        return self._dark_shot_noise

    @dark_shot_noise.setter
    def dark_shot_noise(self, value):
        if value is None:
            self._dark_shot_noise = None
        elif isinstance(value, detector.ShotNoise):
            self._dark_shot_noise = value
        else:
            raise AttributeError('dark_shot_noise must be a ShotNoise object or None')

    @property
    def dsnu(self):
        """Dark signal fixed-pattern noise (DSNU) model. Multiple models are
        combined when specified as a list. If ``None`` (default), DSNU is not
        applied.

        Returns
        -------
        dsnu : :class:`~lentil.detector.FPN`, list_like or None, optional

        """
        return self._dsnu

    @dsnu.setter
    def dsnu(self, value):
        if value is None:
            self._dsnu = None
        elif isinstance(value, detector.FPN):
            self._dsnu = value
        else:
            raise AttributeError('dsnu must be an FPN object or None')

    @property
    def gain(self):
        """Gain model used to convert electrons to DN

        Returns
        -------
        gain : :class:`~lentil.detector.Gain`

        """
        return self._gain

    @gain.setter
    def gain(self, value):
        if value is None:
            self._gain = None
        elif isinstance(value, detector.Gain):
            self._gain = value
        else:
            raise AttributeError('gain must be a Gain object or None')

    @property
    def read_noise(self):
        """Read noise model applied during the readout of the analog signal. If
        ``None`` (default), read noise is not applied.

        Returns
        -------
        read_noise : :class:`~lentil.detector.ReadNoise` or None, optional

        """
        return self._read_noise

    @read_noise.setter
    def read_noise(self, value):
        if value is None:
            self._read_noise = None
        elif isinstance(value, detector.ReadNoise):
            self._read_noise = value
        else:
            raise AttributeError('read_noise must be a ReadNoise object or None')

    @property
    def offset_fpn(self):
        """Offset fixed pattern noise model. Multiple models are applied
        serially when specified as a list. If ``None`` (default), offset FPN is
        not applied.

        Returns
        -------
        offset_fpn : :class:`~lentil.detector.FPN`, list_like, or None, optional

        """
        return self._offset_fpn

    @offset_fpn.setter
    def offset_fpn(self, value):
        if value is None:
            self._offset_fpn = None
        elif isinstance(value, detector.FPN):
            self._offset_fpn = value
        else:
            raise AttributeError('offset_fpn must be an FPN object or None')

    @property
    def charge_diffusion(self):
        """Charge diffusion model describing how charge leaks into adjacent
        pixels. If ``None`` (default), charge diffusion is not applied.

        Returns
        -------
        charge_diffusion : :class:`~lentil.detector.ChargeDiffusion` or None,
        optional

        """
        return self._charge_diffusion

    @charge_diffusion.setter
    def charge_diffusion(self, value):
        if value is None:
            self._charge_diffusion = None
        elif isinstance(value, detector.ChargeDiffusion):
            self._charge_diffusion = value
        else:
            raise AttributeError('charge_diffusion must be a ChargeDiffusion '
                                 'object or None')

    @property
    def defects(self):
        """Detector defect model. Multiple models are combined when specified as
        a list. If ``None`` (default), defects are not included.

        Returns
        -------
        defects : :class:`~lentil.detector.Defect`, list_like, or None, optional

        """
        return self._defects

    @defects.setter
    def defects(self, value):
        if value is None:
            self._defects = None
        elif isinstance(value, detector.Defect):
            self._defects = value
        elif all(isinstance(defect, detector.Defect) for defect in list(value)):
            self._defects = list(value)
        else:
            raise AttributeError('defects must be a Defect object, a list of '
                                 'Defect objects, or None')

    @property
    def cosmic_rays(self):
        """Cosmic ray model. If ``None`` (default), cosmic rays are not
        included.

        Returns
        -------
        cosmic_rays : :class:`~lentil.cosmic.CosmicRay` or None, optional

        """
        return self._cosmic_rays

    @cosmic_rays.setter
    def cosmic_rays(self, value):
        if value is None:
            self._cosmic_rays = None
        elif isinstance(value, detector.CosmicRay):
            self._cosmic_rays = value
        else:
            raise AttributeError('cosmic_rays must be a CosmicRay object or '
                                 'None')

    def frame(self, flux, ts, wave=None, waveunit='nm', oversample=1,
              pixelate=False, collect_charge=False, window=None,
              warn_saturate=False):
        """Simulate a frame.

        Parameters
        ----------
        flux : array_like
            The flux presented to the sensor. ``flux`` should be an (nwave, nrows,
            ncols) datacube with units of photons/s where the first dimension
            represents wavelength if ``collect_charge`` is True or an
            (nrows, ncols) array with units of electrons/s if
            ``collect_charge`` is False.

        ts : float
            Integration time (seconds)

        wave : array_like, optional
            Wavelengths corresponding to each slice in ``flux``. The length of
            ``wave`` must be equal to the number of samples ``n`` in ``flux``.
            If ``collect_charge = False``, ``wave`` is unused and can be
            ``None``.

        waveunit : str, optional
            Units of ``wave``. Defaults is ``nm``

        oversample : int, optional
            Oversampling factor present in ``flux``. Default is 1.

        pixelate : bool, optional
            If True, ``flux`` is convolved with the pixel MTF before rescaling
            the image to native detector sampling. Default is False.
            :func:`pixelate` should only be used if ``oversample`` >= 2

        collect_charge : bool, optional
            If False (default), the input ``flux`` is assumed to be in
            electrons/s. If True, the input flux is assumed to be a data cube of
            photons/s (where the third dimension is wavelength) so the charge
            collection process (photons -> electrons) is performed as a part of
            frame generation.

        window : array_like or None, optional
            Indices of focal plane array to return given as (r_start, r_end,
            c_start, c_end). This definition follows standard Numpy indexing.

        warn_saturate : bool, optional
            Raise a warning when pixels saturate. Default is False.

        Returns
        -------
        frame : ndarray
            Raw frame

        Notes
        -----
        In cases where dynamic noise (jitter, smear, etc.) is present and needs
        to be applied in an oversampled image space, the following approach to
        frame generation is recommended:

            * Compute the polychromatic image in oversampled space with
              ``flatten = False``
            * Compute the binned detector irradiance at the wavelengths present
              in the polychromatic image stack (see:
              :func:`~lentil.radiometry.Spectrum.bin`)
            * Call :func:`collect_charge` to simulate the charge collection
              process. This method will flatten the image, speeding up the
              application of convolutional dynamic noise models
            * Call :func:`frame` with ``collect_charge = False``. Note that in
              this case, the ``wave`` parameter is unused and can be ``None``.

        """

        # If collect charge is True, we need to make sure the user has supplied
        # a wave and waveunit
        if collect_charge:
            assert wave is not None
            assert waveunit is not None

        # accumulate the light signal
        light_e = self.light(flux, ts, wave, waveunit, oversample, pixelate,
                             collect_charge, window)

        # accumulate the dark signal
        dark_e = self.dark(light_e.shape, ts, window)

        # cosmic rays
        if self.cosmic_rays:
            cosmic_e = self.cosmic_rays(light_e.shape, self.pixelscale, ts)
        else:
            cosmic_e = 0

        # combine the signals
        e = light_e + dark_e + cosmic_e

        # apply charge diffusion
        if self.charge_diffusion:
            e = self.charge_diffusion(e)

        # apply offset fpn
        e = self.offset(e, window)

        # read out the frame
        frame = self.readout(e, window, warn_saturate)

        return frame

    def collect_charge(self, count, wave, waveunit='nm', *args, **kwargs):
        """Convert photon count (or flux) to electron count (or flux) by
        applying the detector's wavelength-dependent quantum efficiency.

        Parameters
        ----------
        count : array_like
            The photons presented to the sensor. Should have shape (nwave,
            nrows, ncols)

        wave : array_like
            Wavelengths corresponding to each slice in ``count``. The length of
            ``wave`` must be equal to the number of samples ``nwave`` in ``count``.

        waveunit : str, optional
            Units of ``wave``. Defaults is ``nm``

        Returns
        -------
        elec : ndarray
            Electron count or flux

        Notes
        -----
        The units of ``count`` don't really matter, as long as the user is aware
        that this method converts photons per whatever to electrons per
        whatever. Whatever is nothing for counts and seconds for flux.

        """
        count = np.asarray(count)
        if count.ndim == 2:
            count = count[np.newaxis, ...]

        qe = self.qe.sample(wave, waveunit=waveunit)
        if qe.ndim == 0:
            qe = qe[..., np.newaxis]
        return np.einsum('ijk,i->jk', count, qe)

    def bias(self, shape=None, window=None):
        """Simulate a bias (zero integration time) frame

        Parameters
        ----------
        shape : array_like or None, optional
            Shape to

        window : array_like or None, optional
            Indices of focal plane array to return given as (r_start, r_end,
            c_start, c_end). This definition follows standard Numpy indexing.

        Returns
        -------
        bias : ndarray

        Notes
        -----
        This method is equivalent to calling :func:`frame` with an integration
        time of 0.

        """
        if shape is None:
            shape = self.shape
        else:
            shape = tuple(np.asarray(shape))

        # we don't collect any charge so we just read out a dark frame with zero
        # integration time
        e = self.dark(shape, 0, window)

        # apply charge diffusion
        if self.charge_diffusion:
            e = self.charge_diffusion(e)

        # apply offset fpn
        e = self.offset(e, window)

        # read out the frame
        frame = self.readout(e, window, warn_saturate=False)

        return frame

    def light(self, flux, ts, wave, waveunit, oversample, pixelate, collect_charge,
              window):

        # integrate
        count = self.integrate(flux, ts)

        # photon to electron conversion - there are two possibilities here:
        # 1. count is a wavelength dependent cube so we call collect_charge to
        #    apply the QE and collapse the cube into a single frame
        # 2. count is already collapsed into a single frame
        if collect_charge:
            light_e = self.collect_charge(count, wave, waveunit, oversample, window)
        else:
            light_e = count

        # pixelate
        if pixelate:
            light_e = self.pixelate(light_e, oversample)
        elif oversample > 1:
            light_e = util.rebin(light_e, oversample)

        # shot noise
        if self.shot_noise:
            light_e = self.shot_noise(light_e)

        # prnu
        if self.prnu:
            light_e = self.prnu(light_e, window)

        # no negative electrons
        light_e[light_e < 0] = 0

        return np.floor(light_e)

    def dark(self, shape, ts, window):

        if ts > 0:
            # dark current
            if self.dark_signal:
                dark_flux = self.dark_signal(shape, window)
            else:
                dark_flux = np.zeros(shape)

            # integrate
            dark_e = self.integrate(dark_flux, ts)

            # apply dsnu
            if self.dsnu:
                if isinstance(self.dsnu, detector.FPN):
                    dark_e = self.dsnu(dark_e, window)
                elif isinstance(self.dsnu, (list, tuple)):
                    for dsnu in self.dsnu:
                        dark_e = dsnu(dark_e, window)
                else:
                    raise TypeError('DSNU must be an FPN object or list-like')

            # shot noise
            if self.dark_shot_noise:
                dark_e = self.dark_shot_noise(dark_e)

        else:
            dark_e = np.zeros(shape)

        # no negative electrons
        dark_e[dark_e < 0] = 0

        return np.floor(dark_e)

    def offset(self, electrons, window):
        # offset fpn
        if self.offset_fpn:
            if isinstance(self.offset_fpn, detector.FPN):
                electrons = self.offset_fpn(electrons, window)
            elif isinstance(self.offset_fpn, (list, tuple)):
                for fpn in self.offset_fpn:
                    electrons = fpn(electrons, window)
            else:
                raise TypeError('offset_fpn must be an OffsetFPN object or list-like')

        return electrons

    def readout(self, electrons, window, warn_saturate):
        # read noise
        if self.read_noise:
            electrons = self.read_noise(electrons)

        # no negative electrons
        electrons[electrons < 0] = 0

        # no partial electrons
        electrons = np.floor(electrons)

        # adc
        frame = self.gain(electrons, window, warn_saturate)

        # include focal plane defects
        if self.defects:
            if isinstance(self.defects, detector.Defect):
                frame = self.defects(frame, window)
            elif isinstance(self.defects, (list, tuple)):
                for defect in self.defects:
                    frame = defect(frame, window)

        return frame

    @staticmethod
    def integrate(img, ts):
        return np.floor(img * ts)

    @staticmethod
    def pixelate(img, oversample):
        pixel_sampling = Pixel()
        img = pixel_sampling(img, oversample)
        return util.rescale(img, 1/oversample, order=3, mode='nearest', unitary=True)


class BayerFPA(FPA):
    """Bayer focal plane array base class.

    """
    def __init__(self, pixelscale, shape, bayer_pattern, **kwargs):
        super().__init__(pixelscale, shape, **kwargs)
        self._bayer_pattern = bayer_pattern

        self._red_qe = kwargs.get('red_qe', None)
        self._green_qe = kwargs.get('green_qe', None)
        self._blue_qe = kwargs.get('blue_qe', None)

    @property
    def red_qe(self):
        """Red channel quantum efficiency

        Returns
        -------
        red_qe : :class:`~lentil.radiometry.Spectrum`

        """
        return self._red_qe

    @property
    def green_qe(self):
        """Green channel quantum efficiency

        Returns
        -------
        green_qe : :class:`~lentil.radiometry.Spectrum`

        """
        return self._green_qe

    @property
    def blue_qe(self):
        """Blue channel quantum efficiency

        Returns
        -------
        blue_qe : :class:`~lentil.radiometry.Spectrum`

        """
        return self._blue_qe

    @property
    def bayer_pattern(self):
        """A description of the detector's Bayer pattern. Each value represents
        the order of the red, green, and blue pixels starting with the upper-
        left corner of the image and moving left-to-right, top-to-bottom.` For
        example, ``'GBRG'`` represents the following Bayer pattern:

        = =
        G B
        R G
        = =

        Returns
        -------
        bayer_pattern : ``str``

        """
        if self._bayer_pattern:
            bayer_pattern = self._bayer_pattern.reshape(4)
            return (str(bayer_pattern[0]) + str(bayer_pattern[1]) +
                    str(bayer_pattern[2]) + str(bayer_pattern[3]))
        else:
            return None

    @bayer_pattern.setter
    def bayer_pattern(self, value):
        self._bayer_pattern = np.array(list(value.upper())).reshape(2, 2)

    @property
    def red_kernel(self):
        """Return red kernel"""
        return np.where(self._bayer_pattern == 'R', 1, 0)

    @property
    def blue_kernel(self):
        """Return blue kernel"""
        return np.where(self._bayer_pattern == 'B', 1, 0)

    @property
    def green_kernel(self):
        """Return green kernel"""
        return np.where(self._bayer_pattern == 'G', 1, 0)

    def collect_charge(self, count, wave, waveunit='nm', oversample=1, window=None):
        """Convert photon count (or flux) to electron count (or flux) by
        applying the detector's wavelength-dependent quantum efficiency.
        Additional processing is performed to apply separate QE curves and
        location masks for the separate red, green, and blue channels.

        Parameters
        ----------
        count : array_like
            The photons presented to the sensor. Should have shape (nwave,
            nrows, ncols)

        wave : array_like
            Wavelengths corresponding to each slice in ``count``. The length of
            ``wave`` must be equal to the first dimension in ``count``.

        waveunit : str, optional
            Units of ``wave``. Defaults is ``nm``

        oversample : int, optional
            Oversampling factor present in ``count``. Default is 1.

        window :  array_like or None, optional
            Indices of focal plane array to use given as (r_start, r_end,
            c_start, c_end). This definition follows standard Numpy indexing.

        Returns
        -------
        elec : ndarray
            Electron count or flux

        Notes
        -----
        * The units of ``count`` don't really matter, as long as the user is
          aware that this method converts photons per <whatever> to electrons
          per <whatever>. <Whatever> is nothing for counts and seconds for flux.
        * It is important to correctly specify :attr:`oversample` so that the
          correct pixel location mask is applied to any oversampled data.

        """
        nrow = count.shape[1]//oversample
        ncol = count.shape[2]//oversample

        # build up the bayer image. we do this one channel at a time, and then
        # accumulate the individual channels into one final frame
        red_qe = self.red_qe.sample(wave, waveunit=waveunit)
        red_mosaic = np.tile(self.red_kernel, (nrow//2, ncol//2))
        red_mosaic = ndimage.zoom(red_mosaic, oversample, order=0, mode='wrap')
        red_e = np.einsum('ijk,i->jk', count, red_qe) * red_mosaic

        green_qe = self.green_qe.sample(wave, waveunit=waveunit)
        green_mosaic = np.tile(self.green_kernel, (nrow//2, ncol//2))
        green_mosaic = ndimage.zoom(green_mosaic, oversample, order=0, mode='wrap')
        green_e = np.einsum('ijk,i->jk', count, green_qe) * green_mosaic

        blue_qe = self.blue_qe.sample(wave, waveunit=waveunit)
        blue_mosaic = np.tile(self.blue_kernel, (nrow//2, ncol//2))
        blue_mosaic = ndimage.zoom(blue_mosaic, oversample, order=0, mode='wrap')
        blue_e = np.einsum('ijk,i->jk', count, blue_qe) * blue_mosaic

        return red_e + green_e + blue_e
