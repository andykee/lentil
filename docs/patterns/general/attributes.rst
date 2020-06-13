System Attributes and Common Data
=================================

Defining system attributes
--------------------------
System attributes should be defined as constants in whichever module makes the most
sense. When other modules require access to an attribute, it should be imported from the
module in which it is defined. In this way, attributes can be used throughout the model
but only need to be changed/updated in one place. For example:

.. code:: python3

    import lentil

    FOCAL_LENGTH = 3
    PM_DIAMETER = 0.3

    class TinyPupil(lentil.Pupil):
        def __init__(self):
            super().__init__(focal_length=FOCAL_LENGTH, diameter=PM_DIAMETER)


Loading Numpy ``.npy`` and ``.npz`` files
-----------------------------------------
Numpy data stored in ``.npy`` and ``.npz`` files should be loaded with the
`numpy.load <https://docs.scipy.org/doc/numpy/reference/generated/numpy.load.html>`_
command.

* ``.npy`` files store single arrays
* ``.npz`` files store multiple arrays in a dictionary-like mapping

We'll update the example above to load in some additional Numpy data:


.. code:: python3

    import numpy as np
    import lentil

    FOCAL_LENGTH = 3
    PM_DIAMETER = 0.3

    AMPLITUDE = np.load('amp.npy')
    with np.load('mask.npz') as data:
        MASK = data['mask']

    pupil = lentil.Pupil(focal_length=FOCAL_LENGTH,
                         diameter=DIAMETER,
                         mask=MASK,
                         amplitude=AMPLITUDE)


Loading ``FITS`` files
----------------------
FITS files are most easily loaded using Astropy. You'll have to figure out which Header
Data Unit (HDU) stores the data you need, but the interface is otherwise fairly
straightforward:

.. code:: python3

    from astropy.io import fits

    hdul = fits.open('opd.fits')
    OPD = hdul[0].data


More details are available in the `astropy.io.fits <http://docs.astropy.org/en/stable/io/fits/index.html>`_
documentation.

Loading MATLAB ``.MAT`` files
-----------------------------
For v6 and v7 to v7.2 MAT files, scipy.io.loadmat does the trick. MATLAB 7.3 format MAT
files are HDF5 and not supported by scipy.io.loadmat but can probably be loaded by some
other HDF5 Python library.

.. code:: python3

    from scipy.io import loadmat

    mat_contents = loadmat('opd.mat')
    OPD = mat_contents('opd')


More details are available in the `scipy.io.loadmat <https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html>`_
documentation.

.. note::

    The easiest way to deal with MAT files is to convert them to Numpy arrays with
    `mat2ndarray <https://github.com/andykee/lentil/blob/master/matlab/mat2ndarray.m>`_,
    get them in to Python, and save them using numpy.save.


Loading ``CSV`` files
---------------------

.. code:: python3

    import numpy as np

    OPD = np.genfromtxt('opd.csv', delimiter=',', ship_header=1)

More details are available in the `numpy.genfromtxt <https://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html>`_
documenation.


