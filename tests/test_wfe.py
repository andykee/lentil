import numpy as np
import lentil


def test_power_spectrum_wfe():
    mask = lentil.util.circlemask((256, 256), 128)
    rms = 50e-9
    phase = lentil.wfe.power_spectrum(mask, pixelscale=1/256, rms=rms,
                                      half_power_freq=5, exp=3)

    assert np.std(phase[np.nonzero(phase)])/rms >= 0.8


def test_translation_defocus():
    mask = lentil.util.circlemask((256, 256), 128)
    f_number = np.random.normal(10)
    translation = np.random.uniform(low=-0.5e-3, high=0.5e-3)

    phase = lentil.wfe.translation_defocus(mask, f_number, translation)
    pv_defocus = translation/(8*f_number**2)

    assert np.isclose(np.abs(pv_defocus), np.abs(np.max(phase) - np.min(phase)))
