
PTYPES = ('none', 'pupil', 'image', 'tilt', 'transform')


def ptype(ptype):
    """Create a plane type object

    Parameters
    ----------
    ptype
        Object to be converted to a plane type object. Can be
        a string, a ptype, or None.
    """
    if ptype is None:
        ptype = 'none'

    if isinstance(ptype, PType):
        return ptype
    else:
        return PType(ptype)


class PType:
    def __init__(self, ptype):
        if ptype not in PTYPES:
            raise TypeError(f"plane type '{ptype}' not understood")
        self._key = ptype

    def __eq__(self, other):
        if self._key == other._key:
            return True
        else:
            return False
        
    def __hash__(self):
        return hash(self._key)
    
    def __repr__(self):
        return f"ptype('{self._key}')"

    def __str__(self):
        return self._key
    