import numpy as np
import lentil

amp = lentil.util.circle((512, 512), 256)
# make sure we're not starting with an already normalized array!
assert np.sum(np.abs(amp)**2) != 1.0


def test_normalize_power():
    amp_norm = lentil.modeltools.normalize_power(amp)
    assert np.isclose(np.sum(np.abs(amp_norm)**2), 1.0)


def test_normalize_power_2():
    amp_norm = lentil.modeltools.normalize_power(amp, power=2)
    assert np.isclose(np.sum(np.abs(amp_norm)**2), 2.0)
