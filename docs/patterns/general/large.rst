Larger Models
=============
A simple Lentil model can be written in a single file but as the model grows,
maintaining all the code in a single file becomes more difficult. It usually makes more
sense to divide the code into several smaller *modules* and interact with the modules
via a single *package*. The model package should be divided into modules in some logical
way. One way to separate the code is:

* A top-level module containing the main model and providing linkage between the
  individual components (primarily the optical model, radiometric model, and detector
  model
* An module that defines the model's optical planes
* A module defining the model's radiometric properties (if this is very simple, it can
  be included in the top-level or planes module)
* A module for the detector. If a model needs to represent multiple detectors, they
  can be grouped in a single ``detector.py`` module, or be broken out into individual
  modules.
* A ``data`` subdirectory for holding static data. This subdirectory can be further
  subdivided as needed for easier organization.
* Top-level subdirectories for documentation and tests
* A top-level ``scripts`` subdirectory can be useful for holding commonly used scripts
  for interacting with the model

A fairly standard directory structure might look something like this:

.. code-block:: none

    ~/dev/tiny-lentil
    ├── tiny_telescope/
    │   ├── __init__.py
    │   ├── detector.py
    │   ├── tiny.py
    │   ├── planes.py
    │   ├── radiometry.py
    │   └── data/
    │       ├── detector_qe.csv
    │       ├── pupil_mask.npy
    │       └── sensitivities.npz
    ├── docs/
    ├── scripts/
    ├── tests/
    ├── .gitignore
    ├── README.md
    └── setup.py

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
        name='tiny-monocle',
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
