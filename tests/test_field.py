import numpy as np
import pytest

import lentil.field
from lentil.field import Field

@pytest.mark.parametrize('a_data, a_offset, b_data, b_offset, data, offset', [
    (1, [0,0], 1, [0,0], 1+0j, [0,0]),
    (1, [0,0], 1, [1,0], 0+0j, [0,0]),
    (np.ones((2,2)), [0,0], np.ones((2,2)), [0,0], np.ones((2,2)), [0,0]),
    (np.ones((2,2)), [0,0], np.ones((2,2)), [4,4], 0+0j, [0,0]),
    (1, [0,0], np.ones((3,3)), [-5,-5], np.ones((3,3)), [-5,-5]),
    (np.ones((5,4)), [-2,-2], np.ones((3,3)), [0,-1], np.ones((2,2)), [0,-1])

])
def test_multiply(a_data, a_offset, b_data, b_offset, data, offset):
    a = Field(data=a_data, pixelscale=1, offset=a_offset)
    b = Field(data=b_data, pixelscale=1, offset=b_offset)
    c = a * b
    assert np.array_equal(c.data, data)
    assert np.array_equal(c.offset, offset)


def test_merge():
    a = Field(data=np.ones((5,4)), pixelscale=1, offset=[-2,-2])
    b = Field(data=np.ones((3,3)), pixelscale=1, offset=[0,-1])
    c = lentil.field.merge(a, b)

    c_data = np.zeros((6,5))
    c_data[0:5,0:4] += 1
    c_data[3:6,2:5] += 1

    assert np.array_equal(c.data, c_data)
    assert np.array_equal(c.offset, [-1, -2])
