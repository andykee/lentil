Simple Models
=============
All but the simplest models will benefit from being developed in Python modules. Simple
models not requiring extensive custom code development can often be written in a single
module.

A simple model may have a directory structure that looks something like this:

.. code-block:: none

    ~/dev/tiny-lentil
    ├── tiny_telescope/
    │   ├── __init__.py
    │   ├── tiny.py
    │   ├── detector_qe.csv
    │   ├── pupil_mask.npy
    │   └── sensitivities.npz
    ├── .gitignore
    ├── README.md
    └── setup.py

* The ``tiny.py`` module should contain all of the code needed to represent model
  behaviors
* Additional static files are included as needed

If the ``tiny.py`` module contains the top level object called ``TinyTelescope``, the
``__init__.py`` should look like this:

.. code-block:: python3

    from tiny_telescope.tiny import TinyTelescope

With a model defined in this way, users can interact with it in a standard way:

.. code-block:: python3

    >>> from tiny_telescope import TinyTelescope
    >>> t = TinyTelescope()

**setyp.py**

.. code-block:: python3

    import os
    from setuptools import setup

    setup(
        name='tiny-lentil',
        version='1.0.0',
        author='Tim Apple',
        author_email='tim@apple.com',
        packages=['tiny_telescope'],
        package_data={'tiny_telescope': ['data/*']},
        install_requires=['lentil>=0.1'],
        python_requires='>=3.7',
        )

**.gitignore**

.. code-block:: none

    __pycache__/
    *.pyc
    *.m~
    .DS_Store
    *.swp
