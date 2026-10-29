Contributing
============

This section covers the simplest way to get started working with Lentil
locally so that you can contribute.

Workflow for contributing changes
---------------------------------
We follow a `typical GitHub workflow
<https://docs.github.com/en/get-started/using-github/github-flow>`_ of:

* create a personal fork and local clone of this repo
* create a branch for your contribution
* open a pull request
* work through code review

The sections below cover the contribution steps in more detail.

Create a new development environment
------------------------------------
This is optional, but it's best to start with a fresh development environment
so that you've isolated the packages that you're using for this repository.

To do so, use a tool like `conda <https://docs.conda.io/en/latest/>`_, 
`mamba <https://mamba.readthedocs.io/en/latest/>`_, or 
`venv <https://docs.python.org/3/library/venv.html>`_.

Clone the repository locally
----------------------------
First off you'll need a copy of the ``lentil`` source code. You can clone it
for local development like so:

1. **Fork the repository** so you have your own copy on GitHub. See 
`the GitHub forking guide <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ for more information.
2. **Clone the repository locally** so that you have a local copy to work from:

.. code:: bash
    
    git clone https://github.com/{{ YOUR USERNAME }}/lentil
    cd lentil

Install Lentil from source
--------------------------
Next, install Lentil locally:

.. code:: bash

    pip install -e .

Run the tests
-------------
To run the tests ensure Lentil is instlled, then install `pytest 
<https://docs.pytest.org/en/stable/>`_` with:

.. code:: bash

    pip install pytest

The test suite is run with:

.. code:: bash 

    pytest tests

Build the documentation
-----------------------
To build the documentation, first ensure Lentil is installed, then install
the documentation dependencies with:

.. code:: bash

    cd docs
    pip install -r requirements.txt

The documentation is built with:

.. code:: bash

    make clean && make html

Making releases
---------------
New releases are made whenever significant code changes have been made. Lentil
uses GitHub tags and Actions to automatically publish new releases to PyPi.
The general process for making a new release is:

Prepare the codebase for a new version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Bump ``__version__`` in ``__init__.py``. See the note below about version 
  numbering.
- Make a release commit: ``git commit -m "v0.2.0"``
- Push the commit: ``git push origin main``

.. note::

    Lentil uses `semantic versioning <https://semver.org>`_ to decide whether
    a release is a major, minor, or patch bump. Until Lentil reaches ``v1.0``,
    minor versions are assumed to be breaking and patch versions are assumed
    to include features and/or patches.

Make the release
~~~~~~~~~~~~~~~~
- `Start a new GitHub release <https://github.com/andykee/lentil/releases/new>`_
    - Call the release the current version, e.g. ``v0.2.0``
    - In the **Choose a Tag:** dropdown, type the release name (e.g., 
      ``v0.2.0``) and click "Create new tag"
    - In the **Target:** dropdown, pin it to the release commit that you've
      just pushed.
    - Generate the automatic release notes.
- Confirm the release process completed
    - The `pypi GitHub Action job
      <https://github.com/andykee/lentil/blob/master/.github/workflows/pypi.yml>`_
      completed successfully in the `actions tab 
      <https://github.com/andykee/lentil/actions>`_.
    - `The PyPI version is updated <https://pypi.org/project/lentil>`_


