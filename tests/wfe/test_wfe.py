import numpy as np
import lentil


def test_translation_defocus():
    mask = lentil.util.circlemask((256, 256), 128)
    f_number = np.random.normal(10)
    translation = np.random.uniform(low=-0.5e-3, high=0.5e-3)

    phase = lentil.wfe.translation_defocus(mask, f_number, translation)
    pv_defocus = translation/(8*f_number**2)

    assert np.abs(pv_defocus) == np.abs(np.max(phase) - np.min(phase))
