import numpy as np

import lentil


def test_blackbody_planck_radiance():
    b = lentil.radiometry.Blackbody(wave=np.arange(500, 5000), temp=3000, waveunit='nm',
                                    valueunit='photlam')

    planck_radiance = lentil.radiometry.planck_radiance(wave=np.arange(500, 5000),
                                                        temp=3000, waveunit='nm',
                                                        valueunit='photlam')

    assert np.array_equal(b.value, planck_radiance)


def test_blackbody_add():
    b = lentil.radiometry.Blackbody(wave=np.arange(500, 5000), temp=3000, waveunit='nm',
                                    valueunit='wlam')
    r = b + 10
    assert isinstance(r, lentil.radiometry.Spectrum)
    assert not isinstance(r, lentil.radiometry.Blackbody)
    assert np.array_equal(r.value, b.value + 10)


def test_blackbody_mul():
    b = lentil.radiometry.Blackbody(wave=np.arange(500, 5000), temp=3000, waveunit='nm',
                                    valueunit='photlam')
    r = b * 3
    assert isinstance(r, lentil.radiometry.Spectrum)
    assert not isinstance(r, lentil.radiometry.Blackbody)
    assert np.array_equal(r.value, b.value * 3)
