
class Cache:
    """A dict-like object for caching data."""

    def __init__(self):
        self._cache = {}

    def set(self, key, value):
        """Set a cache value.

        Parameters
        ----------
        key : string

        value : object

        """
        self._cache[key] = value

    def get(self, key, default=None):
        """Get a cache value. If the object doesn't exist in the cache, the
        default value is returned.

        Parameters
        ----------
        key : string

        default : object, optional

        Returns
        -------
        value
            Cached value

        Note
        ----
        The literal value `None` should not be stored in a cache because it is
        impossible to distinguish between the stored value `None` and a cache
        miss signified by a returned value of `None`.

        """
        try:
            return self._cache[key]
        except KeyError:
            return default

    def add(self, key, value):
        """Set a cache value only if it doesn't already exist in the cache.

        Parameters
        ----------
        key : string

        value

        """
        if key not in self.keys():
            self._cache[key] = value

    def delete(self, key):
        """Delete a cache value

        Parameters
        ----------
        key : string

        """
        try:
            self._cache.pop(key)
        except KeyError:
            pass

    def clear(self):
        """Delete all the cache values."""
        self._cache = {}

    def keys(self):
        return self._cache.keys()
