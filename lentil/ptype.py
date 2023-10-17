
PTYPES = ('none', 'pupil', 'image', 'transform')


def ptype(ptype):
    """Create a plane type object

    Parameters
    ----------
    ptype
        Object to be converted to a plane type object.
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
        self._type = ptype

    def __eq__(self, other):
        if self._type == other._type:
            return True
        else:
            return False
        
    def __repr__(self):
        return f"ptype('{self._type}')"

    def __str__(self):
        return self._type
    