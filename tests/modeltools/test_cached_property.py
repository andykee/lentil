import lentil


def test_cached_property_no_cache():
    p = lentil.Pupil(diameter=1, focal_length=1, pixelscale=1, phase=5)
    assert p.phase == 5


def test_cached_property_cache():
    p = lentil.Pupil(diameter=1, focal_length=1, pixelscale=1, phase=5)
    p.cache.set('phase', 12)
    assert p.phase == 12
