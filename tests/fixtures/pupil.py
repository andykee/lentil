import pytest

import lentil

@pytest.fixture
def pupil():

    def _pupil(focal_length, diameter, shape, radius, coeffs=None):

        amplitude = lentil.circle(shape, radius)
        pixelscale = diameter/(2*radius)

        if coeffs is not None:
            opd = lentil.zernike_compose(amplitude, coeffs)
        else:
            opd = 0

        p = lentil.Pupil(amplitude=amplitude,
                         phase=opd,
                         pixelscale=pixelscale,
                         focal_length=focal_length)
        
        return p
    
    return _pupil