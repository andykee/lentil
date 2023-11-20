.. _development.install:

******************************
Install Lentil for development
******************************


Forking the Lentil repo
=======================
Matplotlib is hosted at `andykee/lentil.git <https://github.com/andykee/lentil>`_. 
If you plan on solving issues or submitting pull requests to the main Lentil 
repository, you should first fork this repository by visiting 
`andykee/lentil.git <https://github.com/andykee/lentil>`_ and clicking on the 
``Fork``  button on the top right of the page. See the 
`GitHub documentation <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_
for more details.

Installing from source
======================
If you want to build from source in order to work on Lentil itself, first 
clone the Lentil repository:

.. code-block:: bash

    git clone https://github.com/andykee/lentil.git

If you forked the Lentil repository, you should clone the Lentil repository 
from your fork instead (replacing ``<your-username>`` with your GitHub 
username):

.. code-block:: bash

    git clone https://github.com/<your-username>/lentil.git

Now you can install Lentil in editable mode from the ``lentil`` directory:

.. code-block:: bash

    pip install -e .

Development dependencies
========================
Lentil uses the `pytest <https://doc.pytest.org/en/latest/>`_ framework for
testing. Install it with

.. code-block:: bash

    pip install pytest

The additional Python packages required to build the documentation are
listed in ``docs/requirements.txt`` and can be installed using

.. code-block:: bash

    pip install -r docs/requirements.txt

Running the test suite
======================
To run the tests, in the root directory of your development repository run:

.. code-block:: bash

    pytest tests


Building the documentation
==========================
The documentation source is found in the ``docs/`` directory. The 
configuration file for Sphinx is ``docs/conf.py``. It controls which 
directories Sphinx parses, how the docs are built, and how the extensions are 
used. To build the documentation in html format, cd into ``docs/`` and run:

.. code-block:: bash

    make html

The built docs will be placed in the folder ``docs/_build/html``.