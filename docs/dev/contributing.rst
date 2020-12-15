.. _contributing:

**********************
Contributing to Lentil
**********************

.. note::

    This contribution guide is very heavily influenced by (and in many places
    directly copied from) the `Pandas development guide
    <https://pandas.pydata.org/docs/development/index.html>`_.

Bug reports and enhancement requests

Bug reports and enhancement requests should be filed using Lentil's
`issue tracker <https://github.com/andykee/lentil/issues>`__.

Working with the source code
============================

Version control, Git, and GitHub
--------------------------------
Lentil's source code is hosted on `GitHub <https://github.com/andykee/lentil>`_.
To contribute you'll need `an account <https://github.com/signup/free>`_. Local
version control is handled by `Git <https://git-scm.com/>`_.

`GitHub has instructions <https://help.github.com/set-up-git-redirect>`_ for
installing Git, setting up your SSH key, and configuring Git. All these steps
need to be completed before you can work seamlessly between your local repository
and GitHub.

Some useful resources for learning Git:

* `Github Docs <https://docs.github.com/en>`_
* Matthew Brett's `Pydagouge <http://matthew-brett.github.io/pydagogue/git.html>`_
  notes on Git
* `Oh Shit, Git? <https://ohshitgit.com>`_ for when things go horribly wrong

.. _contributing.forking:

Forking
-------
You will need your own fork to work on the code. Go to the Lentil GitHub page and
hit the ``Fork`` button. You will want to clone your fork to your machine::

    git clone https://github.com/your-user-name/lentil.git lentil-yourname
    cd lentil-yourname
    git remote add upstream https://github.com/andykee/lentil.git

This creates the directory `lentil-yourname` and connects your repository to the
upstream (main project) Lentil repository.

Creating a Python environment
-----------------------------
To test out code changes, you’ll need to build install from source, which
requires a suitable Python environment. To create an isolated Lentil development
environment:

* Install either `Anaconda <https://www.anaconda.com/download/>`_ or `miniconda
  <https://conda.io/miniconda.html>`_
* Make sure your conda is up to date (``conda update conda``)
* Make sure that you have :ref:`cloned the repository <contributing.forking>`
* ``cd`` to the Lentil source directory

We can now create a development environment and install Lentil::

    # Create and activate the build environment:
    conda env create -f environment.yml
    conda activate lentil-dev

    # or with older versions of Anaconda:
    source activate lentil-dev

    # Install Lentil and its dependencies
    python -m pip install -e . --no-build-isolation --no-use-pep517

You should now be able to import Lentil in your development environment::

    $ python
    >>> import lentil

To view your environments::

    conda info -e

To return to your root environment::

    conda deactivate

See the full conda docs `here <https://conda.pydata.org/docs>`_.

Contributing to the code base
=============================

Creating a branch
-----------------
You want your master branch to reflect only production-ready code, so create a feature
branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch. Keep any changes
in this branch specific to one bug or feature so it is clear what the branch brings to
Lentil. You can have many shiny-new-features and switch in between them using the git
checkout command.

When creating this branch, make sure your master branch is up to date with the latest
upstream master version. To update your local master branch, you can do::

    git checkout master
    git pull upstream master --ff-only

When you want to update the feature branch with changes in master after you created the
branch, check the section on :ref:`updating a PR <contributing.update-pr>`.

.. _contributing.commit-code:

Committing your code
--------------------
Once you've made changes, you can see them by typing::

    git status

If you have created a new file, it is not being tracked by git. Add it by typing::

    git add path/to/file-to-be-added.py

Doing 'git status' again should give something like::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

Finally, commit your changes to your local repository with an explanatory message.
Lentil uses a convention for commit message prefixes and layout.  Here are
some common prefixes along with general guidelines for when to use them:

* ENH: Enhancement, new functionality
* BUG: Bug fix
* DOC: Additions/updates to documentation
* TEST: Additions/updates to tests
* PERF: Performance improvement
* CLN: Code cleanup

The following defines how a commit message should be structured.  Please reference the
relevant GitHub issues in your commit message using #1234.

* a subject line with ``< 80`` chars.
* One blank line.
* Optionally, a commit message body.

Now you can commit your changes in your local repository::

    git commit -m

.. _contributing.push-code:

Squashing commits
-----------------
It's possible to combine (or squash) a number of smaller commits into one larger
commit. This helps to keep the project history more concise and readable. The
easiest wat to squash commits is by using interactive rebase. To consider the
most recent ``n`` commits::

    git rebase -i HEAD~<n>

To instead consider all commits including and after a specific commit::

    git rebase -i <after-this-commit-sha1>

The interactive rebase interface provides additional syntax details.

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits::

    git push origin shiny-new-feature

Here ``origin`` is the default name given to your remote repository on GitHub.
You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above you will see something
like::

    origin  git@github.com:yourname/lentil.git (fetch)
    origin  git@github.com:yourname/lentil.git (push)
    upstream        git://github.com/andykee/lentil.git (fetch)
    upstream        git://github.com/andykee/lentil.git (push)

Now your code is on GitHub, but it is not yet a part of the Lentil project. For that to
happen, a pull request needs to be submitted on GitHub.

Review your code
----------------

When you're ready to ask for a code review, file a pull request. Before you do, once
again make sure that you have followed all the guidelines outlined in this document
regarding code style, tests, and documentation. You should also double check your
branch changes against the branch it was based on:

#. Navigate to your repository on GitHub -- https://github.com/your-user-name/lentil
#. Click on ``Branches``
#. Click on the ``Compare`` button for your feature branch
#. Select the ``base`` and ``compare`` branches, if necessary. This will be ``master`` and
   ``shiny-new-feature``, respectively.

Make a pull request
-------------------

If everything looks good, you are ready to make a pull request.  A pull request is how
code from a local repository becomes available to the GitHub community and can be looked
at and eventually merged into the master version.  This pull request and its associated
changes will eventually be committed to the master branch and available in the next
release.  To submit a pull request:

#. Navigate to your repository on GitHub
#. Click on the ``Pull Request`` button
#. You can then click on ``Commits`` and ``Files Changed`` to make sure everything looks
   okay one last time
#. Write a description of your changes in the ``Preview Discussion`` tab
#. Click ``Send Pull Request``.

This request then goes to the repository maintainers, and they will review
the code.

.. _contributing.update-pr:

Upadting a pull request
-----------------------

Based on the review you get on your pull request, you will probably need to make
some changes to the code. In that case, you can make them in your branch,
add a new commit to that branch, push it to GitHub, and the pull request will be
automatically updated.  Pushing them to GitHub again is done by::

    git push origin shiny-new-feature

Another reason you might need to update your pull request is to solve conflicts
with changes that have been merged into the master branch since you opened your
pull request.

To do this, you need to "merge upstream master" in your branch::

    git checkout shiny-new-feature
    git fetch upstream
    git merge upstream/master

If there are no conflicts (or they could be fixed automatically), a file with a
default commit message will open, and you can simply save and quit this file.

If there are merge conflicts, you need to solve those conflicts. See for
example at https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/
for an explanation on how to do this.
Once the conflicts are merged and the files where the conflicts were solved are
added, you can run ``git commit`` to save those fixes.

If you have uncommitted changes at the moment you want to update the branch with
master, you will need to ``stash`` them prior to updating (see the
`stash docs <https://git-scm.com/book/en/v2/Git-Tools-Stashing-and-Cleaning>`__).
This will effectively store your changes and they can be reapplied after updating.

After the feature branch has been update locally, you can now update your pull
request by pushing to the branch on GitHub::

    git push origin shiny-new-feature

Delete your merged branch (optional)
------------------------------------

Once your feature branch is accepted into upstream, you'll probably want to get rid of
the branch. First, merge upstream master into your branch so git knows it is safe to
delete your branch::

    git fetch upstream
    git checkout master
    git merge upstream/master

Then you can do::

    git branch -d shiny-new-feature

Make sure you use a lower-case ``-d``, or else git won't warn you if your feature
branch has not actually been merged.

The branch will still exist on GitHub, so to delete it there do::

    git push origin --delete shiny-new-feature

.. _run_tests:

Running the test suite
----------------------

.. note::

  Running the tests requires `pytest <https://docs.pytest.org/en/latest/>`_.

The tests can then be run directly inside your Git clone by typing::

    pytest tests

Documenting your code
---------------------
Changes should be reflected in the release notes located in ``CHANGES.rst``. This
file contains an ongoing change log for each release. Add an entry to this file to
document your fix, enhancement or (unavoidable) breaking change. Make sure to include
the GitHub issue number when adding your entry (using ``:issue:`1234``` where ``1234``
is the issue/pull request number).

If your code is an enhancement, it is most likely necessary to add usage examples to
the existing documentation. Further, to let users know when this feature was added,
the ``versionadded`` directive is used. The sphinx syntax for that is::

  .. versionadded:: 1.1.0

This will put the text New in version 1.1.0 wherever you put the sphinx directive. This
should also be put in the docstring when adding a new function or method or a new keyword
argument.

Building the documentation
--------------------------

.. note::

  Building the documentation requires `Sphinx <https://www.sphinx-doc.org/en/master/>`_
  and the `PyData Sphinx Theme
  <https://pydata-sphinx-theme.readthedocs.io/en/latest/index.html>`_.

To build the documentation, navigate to your local ``docs/`` directory and run::

    make html

The HTML documentation will be written to ``docs/_build/html``.

If you want to do a full clean build, do::

    make clean && make html

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
