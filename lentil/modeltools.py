import numpy as np

__all__ = ['cached_property', 'normalize_power', 'iterable_amplitude', 'iterable_mask',
           'iterable_phase', 'iterable_segmask']


class cached_property:
    """Property that checks an object's cache for a cached value stored with the
    same attribute name before returning whatever value the in-place property
    defines.

    Note
    ----
    :class:`cached_property` is a drop-in replacement for Python's ``property``
    decorator. It has no side-effects and is safe to use even if Monocle's
    caching backend is not used.

    Example
    -------
    First we'll define a simple example class with a ``surface`` property
    using ``@cached_property``:

    .. code:: python3

        import lentil as le

        class CustomPlane(le.Plane):
            def __init__(self, surface):
                self._surface = surface

            @cached_property
            def surface(self):
                    return 10e-9

    Now we'll take a look at how a cached property behaves:

    .. code:: pycon

        >>> p = CustomPlane()
        >>> p.surface
        1e-08
        >>> p.cache.set('surface', 5e-6)
        >>> p.surface
        5e-06
        >>> p.cache.delete('surface')
        >>> p.surface
        1e-08

    """

    def __init__(self, fget=None, fset=None, fdel=None):
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.__doc__ = getattr(fget, '__doc__')

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, cls=None):
        if instance.cache.get(self.name) is not None:
            return instance.cache.get(self.name)
        else:
            return self.fget(instance)

    def __set__(self, instance, value):
        if self.fset is None:
            raise AttributeError("can't set attribute")
        self.fset(instance, value)

    def setter(self, fset):
        return type(self)(self.fget, fset, self.fdel)


def normalize_power(array, power=1):
    r"""Normalizie the power in an array.

    The total power in an array is given by

    .. math::

        P = \sum{\left|\mbox{array}\right|^2}

    A normalization coefficient is computed as

    .. math::

        c = \sqrt{\frac{p}{\sum{\left|\mbox{array}\right|^2}}}

    The array returned by a will be scaled by the normalization coefficient so
    that its power is equal to :math:`p`.

    Parameters
    ----------
    array : array_like
        Array to be normalized

    power : float, optional
        Desired power in normalized array. Default is 1.

    Returns
    -------
    array : ndarray
        Normalized array

    """
    array = np.asarray(array)
    return array * np.sqrt(power/np.sum(np.abs(array)**2))


def iterable_amplitude(iterable):
    """Construct a common amplitude array from an iterable of planes

    Parameters
    ----------
    iterable : list_like
        List of planes

    Returns
    -------
    amplitude : ndarray

    """
    amplitude = np.array(1)
    for plane in iterable:
        amplitude = amplitude * plane.amplitude
    return amplitude


def iterable_phase(iterable):
    """Construct a common phase array from an iterable of planes

    Parameters
    ----------
    iterable: list_like
        List of planes

    Returns
    -------
    phase: ndarray

    """

    phase = np.array(0)
    for plane in iterable:
        phase = phase + plane.phase
    return phase


def iterable_mask(iterable):
    """Construct a common mask array from an iterable of planes

    Parameters
    ----------
    iterable: list_like
        List of planes

    Returns
    -------
    mask: ndarray

    """
    mask = np.array(1)
    for plane in iterable:
        mask = mask * plane.mask
    return mask


def iterable_segmask(iterable):
    """Construct a common segmask array from an iterable of planes

    Parameters
    ----------
    iterable: list_like
        List of planes

    Returns
    -------
    segmask: ndarray

    """
    if all(plane.segmask is None for plane in iterable):
        return None
    else:
        segmask = np.ones_like(iterable[0].segmask)
        for plane in iterable:
            segmask = segmask * plane.segmask
        return segmask
