
PTYPES = ('none', 'pupil', 'image', 'transform')

class ptype:
    def __init__(self, ptype):
        if ptype is None:
            ptype = 'none'
        
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