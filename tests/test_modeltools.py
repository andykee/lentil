import numpy as np
import lentil


# normalize_power --------------------------------------------------------------
amp = lentil.util.circle((512, 512), 256)


def test_amp_norm():
    # make sure we're not starting with an already normalized array!
    assert np.sum(np.abs(amp)**2) != 1.0


def test_normalize_power():
    amp_norm = lentil.modeltools.normalize_power(amp)
    assert np.isclose(np.sum(np.abs(amp_norm)**2), 1.0)


def test_normalize_power_2():
    amp_norm = lentil.modeltools.normalize_power(amp, power=2)
    assert np.isclose(np.sum(np.abs(amp_norm)**2), 2.0)


# plane iterables --------------------------------------------------------------
class RandomPlane(lentil.Plane):
    def __init__(self):
        super().__init__()
        self.amplitude = np.random.uniform(low=-0.5, high=0.5, size=(10, 10))
        self.phase = np.random.uniform(low=-0.5, high=0.5, size=(10, 10))
        self.mask = np.random.uniform(low=-0.5, high=0.5, size=(10, 10))
        self.segmask = np.random.uniform(low=-0.5, high=0.5, size=(2, 10, 10))


def test_iterable_amplitude():
    p1 = RandomPlane()
    p2 = RandomPlane()
    iterable_amplitude = lentil.modeltools.iterable_amplitude([p1, p2])
    assert np.array_equal(iterable_amplitude, p1.amplitude*p2.amplitude)


def test_iterable_amplitude_default():
    p1 = lentil.Plane()
    p2 = RandomPlane()
    iterable_amplitude = lentil.modeltools.iterable_amplitude([p1, p2])
    assert np.array_equal(iterable_amplitude, p1.amplitude*p2.amplitude)


def test_iterable_phase():
    p1 = RandomPlane()
    p2 = RandomPlane()
    iterable_phase = lentil.modeltools.iterable_phase([p1, p2])
    assert np.array_equal(iterable_phase, p1.phase+p2.phase)


def test_iterable_phase_default():
    p1 = lentil.Plane()
    p2 = RandomPlane()
    iterable_phase = lentil.modeltools.iterable_phase([p1, p2])
    assert np.array_equal(iterable_phase, p1.phase+p2.phase)


def test_iterable_mask():
    p1 = RandomPlane()
    p2 = RandomPlane()
    iterable_mask = lentil.modeltools.iterable_mask([p1, p2])
    assert np.array_equal(iterable_mask, p1.mask*p2.mask)


def test_iterable_mask_default():
    p1 = lentil.Plane()
    p2 = RandomPlane()
    iterable_mask = lentil.modeltools.iterable_mask([p1, p2])
    assert np.array_equal(iterable_mask, p1.mask*p2.mask)


def test_iterable_segmask():
    p1 = RandomPlane()
    p2 = RandomPlane()
    iterable_segmask = lentil.modeltools.iterable_segmask([p1, p2])
    assert np.array_equal(iterable_segmask, p1.segmask*p2.segmask)


def test_iterable_segmask_default():
    p1 = lentil.Plane()
    p2 = RandomPlane()
    iterable_segmask = lentil.modeltools.iterable_segmask([p1, p2])
    assert np.array_equal(iterable_segmask, p2.segmask)


def test_iterable_segmask_all_none():
    p1 = lentil.Plane()
    p2 = lentil.Plane()
    iterable_segmask = lentil.modeltools.iterable_segmask([p1, p2])
    assert iterable_segmask is None
