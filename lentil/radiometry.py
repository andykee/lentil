import copy

import numpy as np
import scipy

# Charlotte - We love you and will miss you forever. You were the light of our lives.

__all__ = ['Spectrum', 'Blackbody', 'Material', 'path_emission', 'path_transmission',
           'planck_radiance', 'planck_exitance', 'vegaflux']

# Constants
H = 6.62606957e-34  # Planck constant, J/s (kg m^2/s^2)
C = 299792456      # Speed of light, m/s
K = 1.3806488e-23  # Boltzmann constant, J/K (kg m^2/K)


class Spectrum:
    """Class for representing spectral quantities.

    :class:`~lentil.radiometry.Spectrum` provides core functionality for
    defining and interacting with spectral data in lentil. Some basic
    arithmetic functions are overloaded, so you can add, multiply, and
    exponentiate :class:`Spectrum` objects

        * with each other
        * by constants
        * with arrays having the same length as :attr:`value`

    Parameters
    ----------
    wave : array_like
        Array of wavelengths.
    value : array_like
        Array of values corresponding to wavelengths in :attr:`wave`.
    waveunit : str, optional
        Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.
    valueunit : str or None, optional
        Value units, as accepted by :func:`Unit`. Default is None.

    Notes
    -----
    When performing arithmetic operations between two
    :class:`Spectrum` objects, the code will behave differently depending on
    the ``wave`` vectors of each :class:`Spectrum`:

        1. If the :attr:`wave` vectors of each :class:`Spectrum` are the same,
        the arithmetic operation is directly performed element-wise on the
        :attr:`value` vectors of each :class:`Spectrum`.

        2. If the :attr:`wave` vectors of each :class:`Spectrum` are
        different, the :attr:`value` vectors of each :class:`Spectrum` are
        interpolated to the same sampling, and then the arithmetic operation is
        performed element-wise on the resulting interpolated :attr:`value`
        vectors. Before interpolating, the :attr:`wave` and :attr:`value`
        vectors are zero-padded to cover the entire wave range represented by
        the union of both :class:`Spectrum` objects. The sampling is determined
        by the smaller sampling of the two :attr:`wave` vectors.

    You can fine tune the sampling, interpolation, and padding methods by
    interfacing directly with the :func:`Spectrum.add` and
    :func:`Spectrum.multiply` methods.

    """

    # TODO: astropy unit cutover
    #   * manage data and units separately
    #   * in __init__, unpack wave and value arrays from unit if passed as a
    #     quantity
    #   * in the wave/value getter, repack as a quantity
    #   * for speed, underlying computations should operate on _wave or _value
    #     and keep the units in mind

    def __init__(self, wave, value, waveunit='nm', valueunit=None):
        self.wave = wave
        self.value = value
        self.waveunit = waveunit
        self.valueunit = valueunit
        self.__array_priority__ = 1.0

        if self._wave.shape != self._value.shape:
            raise ValueError('Wave and value must have the same shape')

    def __add__(self, other):
        return self.add(other)

    def __mul__(self, other):
        return self.multiply(other)

    def __pow__(self, other):
        return self.power(other)

    def __sub__(self, other):
        return self.subtract(other)

    def __truediv__(self, other):
        return self.divide(other)

    __rmul__ = __mul__

    @property
    def wave(self):
        return self._wave

    @wave.setter
    def wave(self, value):
        value = np.asarray(value)

        if np.any(value <= 0):
            raise ValueError('Wavelength values must be greater than zero')

        if not np.all(np.sort(value) == value):
            raise ValueError('Wavelength values must be monotonically increasing')

        if np.any(value[1:] - value[:-1] == 0):
            raise ValueError('Wavelength values must be unique')

        self._wave = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        value = np.asarray(value)
        self._value = value

    @property
    def waveunit(self):
        if self._waveunit is not None:
            return self._waveunit.name
        else:
            return None

    @waveunit.setter
    def waveunit(self, waveunit):
        self._waveunit = Unit(waveunit)

    @property
    def valueunit(self):
        if self._valueunit is not None:
            return self._valueunit.name
        else:
            return None

    @valueunit.setter
    def valueunit(self, valueunit):
        self._valueunit = Unit(valueunit)

    def _ufunc(self, ufunc, other, sampling='min', method='linear', fill_value=0):

        if isinstance(other, (int, float, list, tuple, np.ndarray)):
            wave = self.wave
            try:
                value = ufunc(self.value, other)
            except ValueError:
                raise

        elif isinstance(other, Spectrum):
            wave, self_value, other_value = _interp_common(self, other, sampling,
                                                           method, fill_value)
            value = ufunc(self_value, other_value)
        else:
            raise TypeError(f"can't {ufunc.__name__} Spectrum with object of type "
                            f"{type(other).__name__}")

        return Spectrum(wave, value, self.waveunit, self.valueunit)

    def add(self, other, sampling='min', method='linear', fill_value=0):
        """Add Spectrum and other, element-wise

        Parameters
        ----------
        other : scalar, sequence, or Spectrum
            Any accepted data structure
        sampling : {'min', 'left', 'right', float}, optional
            Method used for computing wavelength sampling when Spectrum
            resampling is required prior to adding

            * 'min' uses the minimum sampling in either Spectrum (Default)
            * 'left' uses the minimum sampling in the left Spectrum
            * 'right' uses the minimum sampling in the right Spectrum
            * <float> uses the float value specified
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (Default)
            * 'quadratic' uses second-order spline interpolation
            * 'cubic' uses third-order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

        Returns
        -------
        result : Spectrum
            Result of the arithmetic operation.

        """
        return self._ufunc(np.add, other, sampling, method, fill_value)

    def multiply(self, other, sampling='min', method='linear', fill_value=0):
        """Multiply Spectrum and other, element-wise

        Parameters
        ----------
        other : scalar, sequence, or Spectrum
            Any accepted data structure
        sampling : {'min', 'left', 'right'} or float, optional
            Method used for computing wavelength sampling when Spectrum
            resampling is required prior to adding

            * 'min' uses the minimum sampling in either Spectrum (Default)
            * 'left' uses the minimum sampling in the left Spectrum
            * 'right' uses the minimum sampling in the right Spectrum
            * <float> uses the float value specified
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (Default)
            * 'quadratic' uses second-order spline interpolation
            * 'cubic' uses third-order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

        Returns
        -------
        result : Spectrum
            Result of the arithmetic operation.

        """
        return self._ufunc(np.multiply, other, sampling, method, fill_value)

    def power(self, other, sampling='min', method='linear', fill_value=0):
        """Spectrum elements raised to powers from other, element-wise

        Parameters
        ----------
        other : scalar, sequence, or Spectrum
            Any accepted data structure
        sampling : {'min', 'left', 'right'} or float, optional
            Method used for computing wavelength sampling when Spectrum
            resampling is required prior to adding

            * 'min' uses the minimum sampling in either Spectrum (Default)
            * 'left' uses the minimum sampling in the left Spectrum
            * 'right' uses the minimum sampling in the right Spectrum
            * <float> uses the float value specified
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (Default)
            * 'quadratic' uses second-order spline interpolation
            * 'cubic' uses third-order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

        Returns
        -------
        result : Spectrum
            Result of the arithmetic operation.

        """
        return self._ufunc(np.power, other, sampling, method, fill_value)

    def subtract(self, other, sampling='min', method='linear', fill_value=0):
        """Subtract Spectrum and other, element-wise

        Parameters
        ----------
        other : scalar, sequence, or Spectrum
            Any accepted data structure
        sampling : {'min', 'left', 'right'} or float, optional
            Method used for computing wavelength sampling when Spectrum
            resampling is required prior to adding

            * 'min' uses the minimum sampling in either Spectrum (Default)
            * 'left' uses the minimum sampling in the left Spectrum
            * 'right' uses the minimum sampling in the right Spectrum
            * <float> uses the float value specified
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (Default)
            * 'quadratic' uses second-order spline interpolation
            * 'cubic' uses third-order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

        Returns
        -------
        result : Spectrum
            Result of the arithmetic operation.

        """
        return self._ufunc(np.subtract, other, sampling, method, fill_value)

    def divide(self, other, sampling='min', method='linear', fill_value=0):
        """Divide Spectrum and other, element-wise

        Parameters
        ----------
        other : scalar, sequence, or Spectrum
            Any accepted data structure
        sampling : {'min', 'left', 'right'} or float, optional
            Method used for computing wavelength sampling when Spectrum
            resampling is required prior to adding

            * 'min' uses the minimum sampling in either Spectrum (Default)
            * 'left' uses the minimum sampling in the left Spectrum
            * 'right' uses the minimum sampling in the right Spectrum
            * <float> uses the float value specified
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (Default)
            * 'quadratic' uses second-order spline interpolation
            * 'cubic' uses third-order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

        Returns
        -------
        result : Spectrum
            Result of the arithmetic operation.

        """
        return self._ufunc(np.divide, other, sampling, method, fill_value)

    def append(self, other, copy=False):
        """Append Spectrum to the end of caller

        Parameters
        ----------
        other : Spectrum
            Spectrum to be appended with
        copy : bool, optional
            If True, the a new :class:`Spectrum` object is returned. If False
            (default), the append operation occurs in place.

        """
        if not isinstance(other, Spectrum):
            raise ValueError('Spectrum objects can only be appended with other '
                             'Spectrum objects')

        if np.any(other.wave <= self.wave):
            raise ValueError()

        if copy:
            new = self.copy()
            new.wave = np.append(new.wave, other.wave)
            new.value = np.append(new.value, other.value)
            return new
        else:
            self.wave = np.append(self.wave, other.wave)
            self.value = np.append(self.value, other.value)

    def copy(self):
        return copy.deepcopy(self)

    def sample(self, wave, method='linear', fill_value=0, waveunit='nm'):
        """Sample the :class:`Spectrum` at a set of desired wavelengths.

        Parameters
        ----------
        wave : array_like or float
            Wavelength set for sampling.
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (default)
            * 'quadratic' uses second order spline interpolation
            * 'cubic' uses third order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested
              points outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.
        waveunit : str
            Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.

        Returns
        -------
        ndarray
            An array of sampled values.

        """
        if waveunit != self.waveunit:
            self.to(waveunit)

        interp = scipy.interpolate.interp1d(self.wave, self.value, kind=method,
                                            copy=False, bounds_error=False,
                                            fill_value=fill_value)

        return interp(wave)

    def resample(self, wave, method='linear', fill_value=0, waveunit='nm'):
        """Sample the :class:`Spectrum` object at a set of desired wavelengths
        and overwrite the object's :attr:`wave` and :attr:`value` attributes.
        This method can be viewed as a destructive version of
        :func:`Spectrum.sample`

        Parameters
        ----------
        wave : array_like or float
            Wavelength set for sampling.
        method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (default)
            * 'quadratic' uses second order spline interpolation
            * 'cubic' uses third order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

            Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.

        """
        self.value = self.sample(wave, method=method, fill_value=fill_value,
                                 waveunit=waveunit)
        self.wave = wave
        self.waveunit = waveunit

    def bin(self, wave, interp_method='simps', ends='symmetric', preserve_power=True,
            sample_method='linear', fill_value=0, waveunit='nm'):
        r"""Compute a binned representation of the :class:`Spectrum` data at a
        given array of central wavelengths.

        The general algorithm computes bin edges, sampling the :class:`Spectrum`
        values at these edges, and numerically integrating within the edges of
        each bin. For :math:`n` bins (centered at the wavelength values in the
        ``wave`` array), there will be :math:`n+1` edges. The edges at each bin
        centered at :math:`\lambda_k` are given by

        .. math::

            e_{\mbox{left}} = \lambda_k - \frac{\lambda_k - \lambda_{k-1}}{2}

            e_{\mbox{right}} = \lambda_k + \frac{\lambda_{k+1} - \lambda_k}{2}

        Special treatment is required for the first and last edges
        (:math:`\lambda_k =` ``wave[0]`` and :math:`\lambda_k =` ``wave[-1]``)
        since :math:`\lambda_{k-1}` and :math:`\lambda_{k+1}` are undefined,
        respectively.

        Two different first/last bin edge methods are available:

        .. image:: /_static/img/spectrum_bin.png

        ========= =============================================== ==================================================================================
        ends      Left edge for ``wave[0]``                       Right edge for ``wave[-1]``
        ========= =============================================== ==================================================================================
        symmetric :math:`\lambda_1-\frac{\lambda_2-\lambda_1}{2}` :math:`\lambda_{\mbox{end}}+\frac{\lambda_{\mbox{end}}-\lambda_{\mbox{end}-1}}{2}`
        inside    :math:`\lambda_1`                               :math:`\lambda_{\mbox{end}}`
        ========= =============================================== ==================================================================================

        Parameters
        ----------
        wave : array_like
            Central wavelength set for binning. At least two wavelengths must be
            provided.
        interp_method : {'simps', 'trapz'}, optional
            Numerical integration method. 'simps' (default) is better for smooth
            data while 'trapz' is better for linear data.
        ends : {'symmetric', 'inside'}, optional
            Method for handling the first and last bin edge value. Default is
            'symmetric'. See the notes below for details on the available
            methods.
        preserve_power : bool, optional

        sample_method : {'linear', 'quadratic', 'cubic'}, optional
            * 'linear' uses linear interpolation (default)
            * 'quadratic' uses second order spline interpolation
            * 'cubic' uses third order spline interpolation
        fill_value : float or array_like, optional
            * If a float, this value will be used to fill in requested points
              outside of the data range.
            * If a two-element array, the first element is used to fill
              `value_new < value[0]` and the second element is used for
              `value_new > value[-1]`.
            * If not provided, the default is 0.

        Returns
        -------
        value : ndarray
            The binned values, same shape as wave.

        """
        wave = np.asarray(wave)

        if wave.size < 2:
            raise ValueError('Spectrum.bin requires a minimum of two wavelengths.\n'
                             'If this Spectrum must be represented by a single '
                             'wavelength, consider using Spectrum.integrate() instead.')

        if waveunit != self.waveunit:
            self.to(waveunit)

        if interp_method == 'trapz':
            dx = np.diff(wave)/2
            x = wave[0:-1] + dx

            if ends == 'symmetric':
                x = np.concatenate([[wave[0]-dx[0]], x, [wave[-1]+dx[-1]]])
            elif ends == 'inside':
                x = np.concatenate([[wave[0]], x, [wave[-1]]])
            else:
                raise ValueError('Unknown ends ', ends)

            # sample
            f = self.sample(x, method=sample_method, fill_value=fill_value)

            # apply the chained trapezoidal rule
            bins = np.array([])
            for k in range(1, f.size):
                bins = np.append(bins, 0.5*(f[k-1]+f[k])*(x[k]-x[k-1]))

        elif interp_method == 'simps':
            dx = np.diff(wave)/2
            wave_mid = wave[0:-1] + dx

            # interleave the wave and wave_mid arrays
            x = np.empty((wave.size + wave_mid.size,), dtype=wave.dtype)
            x[0::2] = wave
            x[1::2] = wave_mid

            if ends == 'symmetric':
                x = np.concatenate([[wave[0]-dx[0]], x, [wave[-1]+dx[-1]]])
            elif ends == 'inside':
                x = np.insert(x, 1, x[0] + (x[1]-x[0])/2)
                x = np.insert(x, -1, x[-1] + (x[-2]-x[-1])/2)
            else:
                raise ValueError('Unknown ends ', ends)

            # sample
            f = self.sample(x, method=sample_method, fill_value=fill_value)

            # apply the chained simpson's rule
            bins = np.array([])
            for k in range(1, x.size, 2):
                bins = np.append(bins, ((x[k+1]-x[k-1])/6) * (f[k-1]+4*f[k]+f[k+1]))

        else:
            raise ValueError('Unknown method ', interp_method)

        if preserve_power:
            norm_factor = self.integrate(np.min(wave), np.max(wave), method=interp_method)/np.sum(bins)
            bins *= norm_factor

        return bins

    def integrate(self, start=None, end=None, method='simps'):
        """Compute the integrated value between ``start`` and ``end``.

        Parameters
        ----------
        start : float, optional
            Lower wavelength bound. If not specified (default),
            ``min(self.wave)`` is used.
        end : float, optional
            Upper wavelength bound. If not specified (default),
            ``max(self.wave)`` is used.
        method : {'simps', 'trapz'}, optional
            Numerical integration method. 'simps' (default) is better for smooth
            data while 'trapz' is better for linear data.

        Returns
        -------
        result : float
            Integral as approximated by ``interp_method``

        """
        if start is None:
            start = np.min(self.wave)
        if end is None:
            end = np.max(self.wave)

        indices = np.intersect1d(np.where(self.wave >= start),
                                 np.where(self.wave <= end))
        wave = self.wave[indices]
        value = self.value[indices]

        if method == 'simps':
            result = scipy.integrate.simpson(x=wave, y=value)
        elif method == 'trapz':
            result = np.integrate.trapezoid(x=wave, y=value)
        else:
            raise ValueError('Unknown method ', method)

        return result

    def ends(self, tol=1e-4):
        """Locate the indices defining the continuous nonzero portion of the
        :attr:`~radiometry.Spectrum.value`

        The first and last indices where the normalized :attr:`Spectrum.value`
        is greater than :attr:`tol` define the bounds of the retained portion of
        the :class:`Spectrum`.

        Parameters
        ----------
        tol : float, optional
            Relative tolerance used to find ends. Default is 1e-4

        """
        if max(self.value) <= 0:
            raise ValueError('Spectrum.trim is only available when '
                             'max(Spectrum.value) > 0')
        normval = self.value/max(self.value)
        index = np.where(normval > tol)
        index_min = index[0][0]
        index_max = index[0][-1]
        return index_min, index_max

    def trim(self, tol=1e-4):
        """Trim the zero or near-zero ends off the :class:`Spectrum` object.

        The first and last indices where the normalized :attr:`Spectrum.value`
        is greater than :attr:`tol` define the bounds of the retained portion of
        the :class:`Spectrum`.

        Parameters
        ----------
        tol : float, optional
            Relative tolerance used to find ends. Default is 1e-4

        Notes
        -----
        If :attr:`~radiometry.Spectrum.value` is all zeros, no trim operation is
        performed and the Spectrum remains unchanged.

        """
        # If the value array is all zeros, do nothing!
        if not self.value.any():
            return

        index_min, index_max = self.ends(tol)
        self.wave = self.wave[index_min:index_max+1]
        self.value = self.value[index_min:index_max+1]

    def crop(self, min_wave, max_wave):
        """Crop a :class:`Spectrum` object."""

        if min_wave > self.wave[0]:
            indx = np.where(min_wave > self.wave)
            self.wave = np.delete(self.wave, indx)
            self.value = np.delete(self.value, indx)

        if max_wave < self.wave[-1]:
            indx = np.where(max_wave < self.wave)
            self.wave = np.delete(self.wave, indx)
            self.value = np.delete(self.value, indx)

    def pad(self, ends, sampling='min', mode='constant', **kwargs):
        """Pad a Spectrum

        Parameters
        ----------
        ends : sequence or array_like, optional
            Wavelength ends to pad to. If None (default), a single value is
            padded to each end of the Spectrum according to the sampling and
            mode parameters.
        sampling : 'min' or float, optional
            * 'min' (Default)
            * float
        mode : {'constant', 'edge'}, optional
            * 'constant' pads with a constant value (Default)
            * 'edge' pads with the edge values
        values : sequence or int, optional
            Used in 'constant'. The (before, after) values to set the padded
            values for each end. int is a shortcut for before = after. Default
            is 0.

        """
        if mode == 'constant':
            if 'values' in kwargs:
                values = np.asarray(kwargs['values'])
                if values.shape == ():
                    values = np.append(values, values)
            else:
                values = (0, 0)
        elif mode == 'edge':
            values = np.array([self.value[0], self.value[-1]])
        else:
            raise ValueError('Unknown mode', mode)

        dwave = _sampling(self.wave, sampling)

        minwave, maxwave = self.wave.min(), self.wave.max()

        ends = np.asarray(ends)

        nleft = int(np.ceil((minwave - ends[0])/dwave)) + 1
        nright = int(np.ceil((ends[1]-maxwave)/dwave)) + 1

        leftwave = np.linspace(ends[0], minwave, nleft)
        leftwave = np.delete(leftwave, -1)
        leftvalue = values[0] * np.ones(leftwave.shape)

        rightwave = np.linspace(maxwave, ends[1], nright)
        rightwave = np.delete(rightwave, 0)
        rightvalue = values[1] * np.ones(rightwave.shape)

        self.wave = np.hstack((leftwave, self.wave, rightwave))
        self.value = np.hstack((leftvalue, self.value, rightvalue))

    def asarray(self):
        """Return :attr:`wave` and :attr:`value` as an array.

        Returns
        -------
        array : ndarray
            Spectral data as an array with ``array[0] = wave`` and ``array[1] = value``.

        """
        return np.array((self.wave, self.value))

    def to(self, *args):
        """Set new wavelength and/or value units. In addition to supporting
        single ``wave`` or ``value`` unit conversion, this method accepts two
        unit arguments to perform simultaneous ``wave`` and ``value`` unit
        conversions.

        Parameters
        ----------
        unit : str
            Unit to convert to, as accepted by :func:`Unit`.

        """
        for unit in args:

            if unit.lower() in ['m', 'um', 'nm', 'angstrom']:
                if self.valueunit in ['photlam', 'flam', 'wlam']:
                    # if the current valueunit are relative to waveunit, we
                    # need to convert to the requested waveunit and update the
                    # valueunit appropriately
                    self.wave = self.wave * self._waveunit.to(unit)
                    self.value = self.value / self._waveunit.to(unit)
                elif self.valueunit is None:
                    # the current valueunit is None, so we assume it is
                    # independent of waveunit and only convert to the
                    # requested waveunit
                    self.wave = self.wave * self._waveunit.to(unit)
                # Update waveunit
                self.waveunit = unit

            elif unit.lower() in ['photlam', 'flam', 'wlam']:
                if self.valueunit is None:
                    raise TypeError("Can't convert from None valueunit to " + unit)
                else:
                    # All conversions are done with wavelength in terms of meters
                    wave = self.wave * self._waveunit.to('meter')
                    value = self.value / self._waveunit.to('meter')

                    # Do the value conversion
                    self.value = self._valueunit.to(value, unit, wave) / Meter().to(self.waveunit)

                    # Update fluxunit
                    self.valueunit = unit

            else:
                raise ValueError('Unknown unit')

    @classmethod
    def from_csv(cls, filename, waveunit, valueunit, zero_negatives=True,
                 header_rows=1):
        """Create a :class:`Spectrum` from a csv file.

        Parameters
        ----------
        filename : str
            csv file to load
        waveunit : str
            Wavelength units, as accepted by :func:`Unit`.
        valueunit : str
            Value units, as accepted by :func:`Unit`.
        zero_negatives : bool, optional
            If true (default), negative values are replaced with zeros.
        header_rows : int, optional
            Number of header rows in csv file to ignore. Default is 1.

        """

        dat = np.genfromtxt(filename, delimiter=',', skip_header=header_rows)
        if zero_negatives:
            dat[dat[:, 1] < 0, 1] = 0

        return cls(dat[:, 0], dat[:, 1], waveunit, valueunit)


def _sampling(wave, method='min'):
    """Compute wave sampling using to specified method.

    Parameters
    ----------
    wave : array_like
        Single wavelength array or 2D collection of wavelength arrays
    method: {'min', 'left', 'right'} or float
        Method used for computing wavelength sampling when Spectrum resampling
        is required prior to adding

        * 'min' uses the minimum sampling in wave
        * 'left' uses the minimum sampling in the left Spectrum
        * 'right' uses the minimum sampling in the right Spectrum
        * <float> uses the float value specified

    """
    if method == 'min':

        if isinstance(wave, (list, tuple)):
            dwave = np.array([])
            for w in wave:
                dwave = np.append(dwave, np.diff(w).min())

            return dwave.min()

        elif isinstance(wave, np.ndarray):
            return np.diff(wave).min()

    elif method == 'left':
        if len(wave) != 2:
            raise ValueError()
        return _sampling(wave[0], method='min')

    elif method == 'right':
        if len(wave) != 2:
            raise ValueError()
        return _sampling(wave[1], method='min')

    elif np.isscalar(method):
        return method

    else:
        raise ValueError('Unknown sampling method', method)


def _intersect(subset, superset):
    """Return the superset indices where the subset overlaps based on its
    subset.min() and subset.max().

    Both subset and superset are assumed to be monotonically increasing

    """
    return np.where((superset >= subset.min()) & (superset <= subset.max()))


def _interp_common(s1, s2, sampling, method, fill_value):
    """Return a common wavelength basis with appropriately interpolated and
    padded value arrays for s1 and s2.

    Parameters
    ----------
    s1, s2 : Spectrum objects

    sampling : {'min', 'left', 'right'} or float
        Method used for computing wavelength sampling when Spectrum resampling
        is required prior to adding

        * 'min' uses the minimum sampling in either Spectrum
        * 'left' uses the minimum sampling in the left Spectrum
        * 'right' uses the minimum sampling in the right Spectrum
        * <float> uses the float value specified

    method : {'linear', 'quadratic', 'cubic'}
        Interpolation method used for sampling Spectrum values.

        * 'linear' uses linear interpolation
        * 'quadratic' uses second-order spline interpolation
        * 'cubic' uses third-order spline interpolation

    fill_value : float
        Value used to fill in requested points outside of the data range.

    Returns
    -------
    commonwave : ndarray

    s1_value : ndarray

    s2_value : ndarray

    """
    # compute a common wavelength array that spans both spectrum and has the
    # desired sampling
    minwave = min(s1.wave.min(), s2.wave.min())
    maxwave = max(s1.wave.max(), s2.wave.max())

    dwave = _sampling((s1.wave, s2.wave), sampling)

    num = int(np.ceil((maxwave - minwave)/dwave))
    commonwave = np.linspace(minwave, maxwave, num + 1)

    # get the portion of commonwave that corresponds to the two spectrum objects
    s1_index = _intersect(s1.wave, commonwave)
    s1_wave = commonwave[s1_index]

    s2_index = _intersect(s2.wave, commonwave)
    s2_wave = commonwave[s2_index]

    # sample each Spectrum at the requested sampling
    s1_samplevalue = s1.sample(s1_wave, method=method, fill_value=fill_value, waveunit=s1.waveunit)
    s2_samplevalue = s2.sample(s2_wave, method=method, fill_value=fill_value, waveunit=s2.waveunit)

    # create nominal value arrays
    s1_value = fill_value * np.ones(commonwave.shape)
    s2_value = fill_value * np.ones(commonwave.shape)

    # insert the sampled values into the appropriate slots
    s1_value[s1_index] = s1_samplevalue
    s2_value[s2_index] = s2_samplevalue

    return commonwave, s1_value, s2_value


class Blackbody(Spectrum):
    """Class for representing a blackbody emitter.

    Parameters
    ----------
    wave : array_like
        Array of wavelengths or wavenumbers
    temp : float
        Temperature in K
    waveunit : str
        Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.
    valueunit : str
        Flux units, as accepted by :func:`Unit`. Default is ``photlam``.

    Examples
    --------
    Create a Blackbody object with wavelength range 400-4000 nm and a temperature of
    4000K:

    .. plot::
        :include-source:
        :scale: 50

        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> import lentil
        >>> wave = np.arange(400,4000)
        >>> temp = 5000
        >>> src = lentil.radiometry.Blackbody(wave,temp,waveunit='nm')
        >>> plt.plot(src.wave, src.value), plt.grid()
        >>> plt.xlabel('Wavelength [nm]'), plt.ylabel('Flux [photons/sec/m^2/sr]')

    """
    def __init__(self, wave, temp, waveunit='nm', valueunit='photlam'):
        wave = np.asarray(wave)
        flux = planck_radiance(wave, temp, waveunit, valueunit)
        self.temp = temp
        self.sample_fn = planck_radiance

        super().__init__(wave, flux, waveunit, valueunit)

    def sample(self, wave, waveunit='nm', *args, **kwargs):
        # sample_fn will be planck_radiance if the Blackbody was created from its
        # standard constructor. If instead, Blackbody was created from the vegamag
        # classmethod, sample_fn will be planck_exitance
        return self.sample_fn(wave, self.temp, waveunit, valueunit=self.valueunit)

    def sample_vegamag(self, wave, temp, waveunit, valueunit):
        # Get Vega zero point data for requested band
        E0, wave0 = vegaflux(self.band, waveunit)

        # Compute exitance of Vega at equivalent wavelength in the observing
        # band
        M0 = planck_exitance(wave0, temp, waveunit, valueunit)

        # Compute the source exitance over the desired wavelengths
        M = planck_exitance(wave, temp, waveunit, valueunit)

        # Scale the source irradiance by the requested magnitude
        E = E0 * (M/M0)*10**(-0.4*self.mag)

        return E

    # TODO: this should really return a different type. when that is done, remove
    # the sample_fn garbage and move sample_vegamag to the new object's sample
    # method (should actually create a flux method that takes band, mag, temp,
    # etc and call that in both the constructor and the sample method)
    @classmethod
    def vegamag(cls, wave, temp, mag, band, waveunit='nm', valueunit='photlam'):
        r""" Create a :class:`Blackbody` with its irradiance scaled to an
        apparent magnitude.

        Parameters
        ----------
        wave : array_like
            Array of wavelengths or wavenumbers
        temp : float
            Temperature in K
        mag : float
            Apparent magnitude
        band : str
            Observing band closest to the wavelengths in :attr:`wave`. See
            :func:`vegaflux` for supported bands.
        waveunit : str
            Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.
        valueunit : str
            Flux units, as accepted by :func:`Unit`. Default is ``photlam``.

        Examples
        --------
        Create a Blackbody object with wavelength range 400-900 nm, a temperature of
        4000K, and a visible magnitude of 2:

        .. plot::
            :scale: 50
            :include-source:

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> import lentil
            >>> wave = np.arange(400,900)
            >>> temp = 5000
            >>> mag = 2
            >>> band = 'V'  # Chosen for best overlap with wavelength range
            >>> src = lentil.radiometry.Blackbody.vegamag(wave,temp,mag,band,waveunit='nm')
            >>> plt.plot(src.wave, src.value), plt.grid()
            >>> plt.xlabel('Wavelength [nm]'), plt.ylabel('Flux [photons/sec/m^2]')

        Notes
        -----
        This class implements a special formulation of a :class:`Blackbody`
        which computes irradiance for a star given a magnitude :math:`mag` in
        an observing band chosen from :func:`vegaflux`. The general calculation
        performed here is as follows:

        1. Look up zero point irradiance :math:`E_{\mbox{Vega}}` and central
           wavelength :math:`\lambda_{0,\mbox{Vega}}` in the desired observing
           band
        2. Compute exitance :math:`M_0` for the desired star at
           :math:`\lambda_{0,\mbox{Vega}}`
        3. Compute exitance :math:`M` for the desired star over the requested
           wavelength range
        4. Compute apparent magnitude scaling factor
           :math:`f = 10^{-0.4\mbox{mag}}`
        5. Compute star irradiance as
           :math:`E = E_{\mbox{Vega}} \ f \left(\frac{M}{M_0}\right)`

        """
        wave = np.asarray(wave)

        # Get Vega zero point data for requested band
        E0, wave0 = vegaflux(band, waveunit)

        # Compute exitance of Vega at equivalent wavelength in the observing
        # band
        M0 = planck_exitance(wave0, temp, waveunit, valueunit='photlam')

        # Compute the source exitance over the desired wavelengths
        M = planck_exitance(wave, temp, waveunit, valueunit='photlam')

        # Scale the source irradiance by the requested magnitude
        E = E0 * (M/M0)*10**(-0.4*mag)

        # Construct a blackbody object
        self = cls(wave, temp, waveunit, valueunit)
        self.value = E

        # Populate this object with additional attributes
        self.mag = mag
        self.band = band
        self.sample_fn = self.sample_vegamag

        return self


class Material:
    """Class for representing a material with radiometric properties.

    Parameters
    ----------
    transmission : :class:`~lentil.radiometry.Spectrum` or float, optional
        Spectral transmission. If not specified, a default transmission is
        created which has no radiometric effect.
    emission : :class:`~lentil.radiometry.Spectrum` or float, optional
        Spectral thermal emission of the element. If not specified, a default
        emission is created which has no radiometric contribution.
    contam : float or :class:`~lentil.radiometry.Spectrum`, optional
        Material contamination factor. Can be used for representing EOL/BOL
        optical properties. Should be within [0,1]. Both
        :attr:`~lentil.radiometry.Material.transmission` and
        :attr:`~lentil.radiometry.Material.emission` are multiplied by contam
        before being returned.

    Notes
    -----
    There is no explicit reflection attribute. Reflective optics should use the
    :attr:`~lentil.radiometry.Material.transmission` attribute to represent
    reflectivity.

    See Also
    --------
    radiometry.path_transmission
    radiometry.path_emission

    """
    def __init__(self, transmission=1, emission=0, contam=1):
        self.transmission = transmission
        self.emission = emission
        self.contam = contam

    @property
    def transmission(self):
        return self.contam * self._transmission

    @transmission.setter
    def transmission(self, value):
        self._transmission = value

    @property
    def emission(self):
        # TODO: make sure this is correct
        return self.contam * self._emission

    @emission.setter
    def emission(self, value):
        self._emission = value


def path_transmission(iterable):
    """Construct a common transmission from an iterable of :class:`Spectrum` or
    :class:`Material` objects or numeric types

    Parameters
    ----------
    iterable: list_like
        List of :class:`Spectrum` or :class:`Material` objects or numeric types

    Returns
    -------
    transmission: :class:`~lentil.radiometry.Spectrum`

    """
    transmission = 1
    for item in iterable:
        if isinstance(item, Material):
            transmission = transmission * item.transmission
        else:
            transmission = transmission * item
    return transmission


def path_emission(iterable, emission=0):
    """Construct a common emission from an iterable of :class:`Material` objects

    Parameters
    ----------
    iterable: list_like
        List of :class:`Material` objects

    emission: scalar or :class:`~lentil.radiometry.Spectrum`, optional
        Emission collected upstream of the iterable path. This could be used to
        represent things like zodiacal background or warm cavities emitting
        into the first element of a path.

    Returns
    -------
    emission: :class:`~lentil.radiometry.Spectrum`

    """
    for item in iterable:
        emission = (emission * item.transmission) + item.emission
    return emission


def planck_exitance(wave, temp, waveunit='nm', valueunit='wlam'):
    r"""Compute the Planck law spectral exitance from a blackbody radiator at
    the given temperature.

    .. math::

        M_{\lambda}(T) = \frac{2\pi hc^2}{\lambda^5 \exp \left(\frac{hc}{\lambda k T}\right)-1}

    Parameters
    ----------
    wave : array_like or float
        Wavelength or array of wavelengths
    temp : float
        Blackbody temperature in K
    waveunit : str
        Wavelength units, as accepted by :func:`Unit`.
    valueunit : str
        Flux units, as accepted by :func:`Unit`.

    Returns
    -------
    ndarray
        Spectral exitance in ``valueunit``.

    """
    # convert wave to meters
    wave = wave * Unit(waveunit).to('meter')

    # compute flux in W m^-2 sr^-1 m^-1
    flux = 2*np.pi*H*C**2/(wave**5*(np.exp(H*C/(wave*K*temp))-1))

    # do flux conversion (if necessary)
    if valueunit == 'wlam':
        # just convert flux back to W m^-2 sr^-1 <waveunit>^-1
        return flux / Meter().to(waveunit)
    else:
        # convert to desired fluxunit and back to waveunit
        return Wlam().to(flux, valueunit, wave) / Meter().to(waveunit)


def planck_radiance(wave, temp, waveunit='nm', valueunit='wlam'):
    r"""Compute the Planck law spectral radiance from a blackbody radiator at
    the given temperature.

    .. math::

        M_{\lambda}(T) = \frac{2hc^2}{\lambda^5 \exp \left(\frac{hc}{\lambda k T}\right)-1}

    Parameters
    ----------
    wave : array_like or float
        Wavelength or array of wavelengths
    temp : float
        Blackbody temperature in K
    waveunit : str
        Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.
    valueunit : str
        Flux units, as accepted by :func:`Unit`. Default is ``wlam``.

    Returns
    -------
    ndarray
        Spectral radiance in ``valueunit sr^-1``.

    """
    # convert wave to meters
    wave = wave * Unit(waveunit).to('meter')

    # compute flux in W m^-2 m^-1
    flux = 2*H*C**2/(wave**5*(np.exp(H*C/(wave*K*temp))-1))

    # do flux conversion (if necessary)
    if valueunit == 'wlam':
        # just convert flux back to W m^-2 <waveunit>^-1
        return flux / Meter().to(waveunit)
    else:
        # convert to desired valueunit and back to waveunit
        return Wlam().to(flux, valueunit, wave) / Meter().to(waveunit)


def vegaflux(band, waveunit='nm', valueunit='photlam'):
    """Get the flux of Vega in a desired observing band.

    ====== ================== =======
    Band   Central Wavelength Flux
    ====== ================== =======
    ``U``  360 nm             1790 Jy
    ``B``  438 nm             4036 Jy
    ``V``  545 nm             3636 Jy
    ``R``  641 nm             3064 Jy
    ``I``  798 nm             2416 Jy
    ``J``  1220 nm            1589 Jy
    ``H``  1630 nm            1021 Jy
    ``K``  2190 nm            640 Jy
    ``W1`` 3353 nm            310 Jy
    ``W2`` 4603 nm            172 Jy
    ``W3`` 11561 nm           31.7 Jy
    ``W4`` 22088 nm           8.36 Jy
    ====== ================== =======

    Parameters
    ----------
    band : str
        Observing band.
    waveunit : str
        Wavelength units, as accepted by :func:`Unit`. Default is ``nm``.
    valueunit : str
        Flux units, as accepted by :func:`Unit`. Default is ``photlam``.

    Returns
    -------
    flux : float
        Vega flux in :attr:`band` expressed in :attr:`valueunit`

    wave : float
        Central wavelength in :attr:`band` expressed in :attr:`waveunit`

    References
    ----------
    [1] Bessell et al. (1998) Johnson-Cousins-Class System

    [2] Wright et al. (2010) Wide-Field Infrared Survey Explorer

    """
    vega = {
        # band    nm                   Jy
        'U':  {'wave': 360e-9,   'flux': 1790},
        'B':  {'wave': 438e-9,   'flux': 4036},
        'V':  {'wave': 545e-9,   'flux': 3636},
        'R':  {'wave': 641e-9,   'flux': 3064},
        'I':  {'wave': 798e-9,   'flux': 2416},
        'J':  {'wave': 1220e-9,  'flux': 1589},
        'H':  {'wave': 1630e-9,  'flux': 1021},
        'K':  {'wave': 2190e-9,  'flux': 640},
        'W1': {'wave': 3353e-9,  'flux': 310},
        'W2': {'wave': 4603e-9,  'flux': 172},
        'W3': {'wave': 11561e-9, 'flux': 31.7},
        'W4': {'wave': 22088e-9, 'flux': 8.36}
        }

    band = band.upper()

    if band not in vega:
        raise ValueError('Unknown band: ' + band)

    wave = vega[band]['wave']
    flux = vega[band]['flux']

    flux = flux * 1e-26                 # Jy         -> W/m^2 Hz
    flux = flux * C/(wave**2)           # W/m^2 Hz   -> W/m^2 m
    flux = flux * wave/(H*C)            # W/m^2 m    -> ph/s m^2 m

    # do flux conversion (if necessary)
    if valueunit == 'photlam':
        # just convert ph/s m^2 m^-1 -> ph/sm^-2 <waveunit>^-1
        flux = flux / Meter().to(waveunit)
    else:
        # convert to desired valueunit and back to waveunit. Note that wave
        # is still in terms of meters here.
        flux = Photlam().to(flux, valueunit, wave) / Meter().to(waveunit)

    wave = wave * Meter().to(waveunit)  # m -> <waveunit>

    return flux, wave


def Unit(name):
    r"""Generate a Unit object.

    **Wavelength Units**

    ===================== ================== ==========================
    Name                  Class Object       Units
    ===================== ================== ==========================
    ``m``, ``meter``      :class:`Meter`     SI base unit
    ``um``, ``micron``    :class:`Micron`    :math:`10^{-6}\ \mbox{m}`
    ``nm``, ``nanometer`` :class:`Nanometer` :math:`10^{-9}\ \mbox{m}`
    ``angstrom``          :class:`Angstrom`  :math:`10^{-10}\ \mbox{m}`
    ===================== ================== ==========================

    **Flux Units**

    =========== ================= ===========================================
    Name        Class Object      Units
    =========== ================= ===========================================
    ``photlam`` :class:`Photlam`  :math:`\mbox{photons s}^{-1} \mbox{m}^{-2}`
    ``wlam``    :class:`Wlam`     :math:`\mbox{W m}^{-2}`
    ``flam``    :class:`Flam`     :math:`\mbox{erg s}^{-1} \mbox{cm}^{-2}`
    =========== ================= ===========================================

    Parameters
    ----------
    name : str
        Wavelength or flux units as accpeted by :func:`Unit`.

    Returns
    -------
    Unit or ``None``
        Unit object. ``None`` means unitless.

    """
    if name is not None:
        if name.lower() in ['meter', 'm']:
            return Meter()
        elif name.lower() in ['um', 'micron']:
            return Micron()
        elif name.lower() in ['nm', 'nanometer']:
            return Nanometer()
        elif name.lower() == 'angstrom':
            return Angstrom()
        elif name.lower() == 'photlam':
            return Photlam()
        elif name.lower() == 'flam':
            return Flam()
        elif name.lower() == 'wlam':
            return Wlam()
        else:
            raise ValueError('Invalid unit: ' + name)
    else:
        return None


class _Unit:
    def __str__(self):
        return self.name


class Meter(_Unit):
    """Meter waveunit."""
    name = 'm'

    @staticmethod
    def to(waveunit):
        """Convert waveunit.

        ===== ============ ======================
        From  To           Conversion
        ===== ============ ======================
        ``m`` ``m``        :math:`1`
        ``m`` ``um``       :math:`10^{-6}` ``m``
        ``m`` ``nm``       :math:`10^{-9}` ``m``
        ``m`` ``angstrom`` :math:`10^{-10}` ``m``
        ===== ============ ======================

        Parameters
        ----------
        waveunit : str
            Wavelength units, as accpeted by :func:`Unit`.

        Returns
        -------
        int
            :attr:`wave` scaling factor

        """
        if waveunit.lower() in ['m', 'meter']:
            return 1
        elif waveunit.lower() in ['um', 'micron']:
            return 1e6
        elif waveunit.lower() in ['nm', 'nanometer']:
            return 1e9
        elif waveunit.lower() == 'angstrom':
            return 1e10
        else:
            raise ValueError('Unknown waveunit: ' + waveunit)


class Micron(_Unit):
    """Micron waveunit."""
    name = 'um'

    @staticmethod
    def to(waveunit):
        """Convert waveunit.

        ====== ============ ======================
        From   To           Conversion
        ====== ============ ======================
        ``um`` ``m``        :math:`10^{-6}` ``um``
        ``um`` ``um``       :math:`1`
        ``um`` ``nm``       :math:`10^{3}` ``um``
        ``um`` ``angstrom`` :math:`10^{4}` ``um``
        ====== ============ ======================

        Parameters
        ----------
        waveunit : str
            Wavelength units, as accepted by :func:`Unit`.

        Returns
        -------
        int
            :attr:`wave` scaling factor

        """
        if waveunit.lower() in ['m', 'meter']:
            return 1e-6
        elif waveunit.lower() in ['um', 'micron']:
            return 1
        elif waveunit.lower() in ['nm', 'nanometer']:
            return 1e3
        elif waveunit.lower() == 'angstrom':
            return 1e4
        else:
            raise ValueError('Unknown waveunit: ' + waveunit)


class Nanometer(_Unit):
    name = 'nm'

    @staticmethod
    def to(waveunit):
        if waveunit.lower() in ['m', 'meter']:
            return 1e-9
        elif waveunit.lower() in ['um', 'micron']:
            return 1e-3
        elif waveunit.lower() in ['nm', 'nanometer']:
            return 1
        elif waveunit.lower() == 'angstrom':
            return 1e1
        else:
            raise ValueError('Unknown waveunit: ' + waveunit)


class Angstrom(_Unit):
    name = 'angstrom'

    @staticmethod
    def to(waveunit):
        if waveunit.lower() in ['m', 'meter']:
            return 1e-10
        elif waveunit.lower() in ['um', 'micron']:
            return 1e-4
        elif waveunit.lower() in ['nm', 'nanometer']:
            return 1e-1
        elif waveunit.lower() == 'angstrom':
            return 1
        else:
            raise ValueError('Unknown waveunit: ' + waveunit)


class Photlam(_Unit):
    name = 'photlam'
    longname = 'photons s^-1 m^-2'

    @staticmethod
    def to(flux, fluxunit, wave):
        if fluxunit.lower() == 'wlam':
            return flux * (H*C)/wave
        elif fluxunit.lower() == 'flam':
            return flux * (H*C)/wave * 1e7 * 1e-4
        elif fluxunit.lower() == 'photlam':
            return flux
        else:
            raise ValueError('Unknown fluxunit: ' + fluxunit)


class Flam(_Unit):
    name = 'flam'
    longname = 'erg s^-1 cm^-2'

    @staticmethod
    def to(flux, fluxunit, wave):
        if fluxunit.lower() == 'photlam':
            return flux * wave/(H * C) * 1e-7 * 1e4
        elif fluxunit.lower() == 'wlam':
            return flux * 1e-7 * 1e4
        elif fluxunit.lower() == 'flam':
            return flux
        else:
            raise ValueError('Unknown fluxunit: ' + fluxunit)


class Wlam(_Unit):
    name = 'wlam'
    longname = 'W m^-2'

    @staticmethod
    def to(flux, fluxunit, wave):
        if fluxunit.lower() == 'photlam':
            return flux * wave/(H*C)
        elif fluxunit.lower() == 'flam':
            return flux * 1e7 * 1e-4
        elif fluxunit.lower() == 'wlam':
            return flux
        else:
            raise ValueError('Unknown fluxunit: ' + fluxunit)
