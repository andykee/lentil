import numpy as np
import lentil.cache


def test_cache_set():
    cache = lentil.cache.Cache()
    value = np.random.random()
    cache.set('value', value)
    assert cache.get('value') == value


def test_cache_delete():
    cache = lentil.cache.Cache()
    value = np.random.random()
    cache.set('value', value)
    cache.delete('value')
    assert 'value' not in cache.keys()


def test_cache_delete_no_key():
    cache = lentil.cache.Cache()
    cache.delete('value')


def test_cache_clear():
    cache = lentil.cache.Cache()
    value = np.random.random()
    cache.set('value', value)
    cache.clear()
    assert len(cache.keys()) == 0
