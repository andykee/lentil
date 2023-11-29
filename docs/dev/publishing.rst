.. _publishing:

Publishing to PyPi
==================

Prepare the release
-------------------
1. Increment the version number appropriately in `lentil/__init__.py
   <https://github.com/andykee/lentil/blob/master/lentil/__init__.py>`_ according to
   `PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_.
2. Document contents of the release in `CHANGES.rst
   <https://github.com/andykee/lentil/blob/master/CHANGES.rst>`_
3. Commit the changes from above on the master branch with a commit message equal to the
   short version name (i.e. ``v0.3.2`` or ``v0.5.0b2``).
4. Push the updated master branch and create a new release.

Build the release
-----------------
.. code:: bash

    $ python setup.py sdist bdist_wheel

Upload the release
------------------
.. code:: bash

    $ twine upload dist/*

.. note::

    Only core-team members are able to publish new releases to PyPi.
