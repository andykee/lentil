.. _user.matlab:

************************
Using Lentil with MATLAB
************************

MATLAB Lentil Interface
=======================

Calling Python libraries from MATLAB is as simple as configuring MATLAB to use the
appropriate Python implementation and prepending the Python command with ``py.``. Before
using Lentil in MATLAB, be sure that your Python version is `compatible with your MATLAB 
release <https://www.mathworks.com/support/requirements/python-compatibility.html>`_ and 
that you have :ref:`configured MATLAB to use the correct version of Python installed on 
your system <user.matlab.config>`. 

A MATLAB helper function is available :download:`here </_static/matlab/check_pyenv.m>`
for verifying compatibility between MATLAB and Python versions.

It's possible to read and write `.npy` files in MATLAB using `npy-matlab
<https://github.com/kwikteam/npy-matlab>`_

Interacting directly with Lentil from MATLAB
--------------------------------------------
It is possible to interact directly with Lentil within MATLAB. For example, we can
create a simple :func:`~circle` mask with:

.. code-block:: matlab

    >> mask = py.lentil.circle([int16(256),int16(256)],int16(128));

Note that MATLAB `automatically does type conversion 
<https://www.mathworks.com/help/matlab/matlab_external/passing-data-to-python.html>`_ 
when MATLAB data is passed to Python but because MATLAB's default numeric type is 
double-precision floating point, the ``shape`` and ``radius`` parameters must
be explicitly cast to ints before passing them to Python.

Developing a MATLAB interface to a Lentil model
-----------------------------------------------
The easiest way to provide a MATLAB interface to an underlying Lentil model is to
develop a separate MATLAB class that mimics the interface of a Lentil class and handles
any data formatting issues that arise between MATLAB and Python. In this way, the model
logic all resides within the Python code, and the MATLAB code is only responsible for
managing the interface between languages. An example following this approach is
available :ref:`here<examples.matlab_interface>`.

Finally, a few links that may be helpful when developing a MATLAB interface:

* `Calling Python functions from MATLAB <https://www.mathworks.com/help/matlab/matlab_external/python-function-arguments.html>`_
* `Passing a Python function kwargs from MATLAB <https://www.mathworks.com/help/matlab/ref/pyargs.html>`_

.. _user.matlab.config:

Configuring MATLAB to use the correct version of Python
=======================================================
MATLAB automatically selects and loads a Python version when you type a Python command.
It commonly defaults to using Python 2.7, which is not at all what we want. We need to
tell MATLAB where to find the correct version of Python (>=3.7). For MATLAB r2020a or
later:

.. code-block:: matlab

    pyenv('Version', '/path/to/python3/executable')

For MATLAB r2019b or earlier:

.. code-block:: matlab

    pyversion '/path/to/python3/executable'

The Python version can be set temporarily by executing the above command when MATLAB
launches or automatically by adding the command to your ``startup.m`` file.

.. note::
    If you're using virtual environments to manage different Lentil models, the
    ``pyenv/pyversion`` configuration specified above won't work. Instead, you'll need
    to call ``pyenv/pyversion`` with the correct virtual environment version before
    working with a model. Don't forget to call the MKL conflict fix as well. This is
    annoying. Sorry about that. It's a MATLAB "feature".

.. warning::
    Once you've set ``pyenv/pyversion`` within a MATLAB session, the only way to change
    it is to restart MATLAB. This means that if you're working with virtual
    environments to manage different models, you'll have to restart MATLAB each time you
    want to switch models. This is annoying. Sorry about that. It's a MATLAB "feature".

For more help on setting MATLAB's Python version, see
`System and Configuration Requirements <https://www.mathworks.com/help/matlab/matlab_external/system-and-configuration-requirements.html>`_.


Resolving MKL Conflicts
-----------------------
MATLAB doesn't always load the correct libraries the underlying Python code relies on.
In particular, there seems to be some confusion about when to load MKL. There is no
telltale sign this has occurred. Sometimes MATLAB crashes while other times Python
method calls will error out with messages that may or may not be useful. The following
command will clear up MATLAB's confusion by handing control of which libraries Python
needs back to Python:

.. code-block:: matlab

    py.sys.setdlopenflags(int32(10));

This command sets the ``RTLD_NOW`` and ``RTLD_DEEPBIND`` flags when the active Python
instance calls ``dlopen()`` `[1]`_ `[2]`_ `[3]`_. Note that this command is Unix only
and must be called before the Python interpreter is loaded within MATLAB but after
``pyenv/pyversion`` is set, making it a prime candidate for inclusion in ``startup.m``.

.. _[1]: https://www.mathworks.com/matlabcentral/answers/327193-calling-python-module-from-matlab-causes-segmentation-fault-in-h5py#answer_296569
.. _[2]: http://man7.org/linux/man-pages/man3/dlopen.3.html
.. _[3]: https://docs.python.org/3.6/library/sys.html#sys.setdlopenflags


Troubleshooting
===============

Debugging MATLAB's Undefined variable "py" or function "py.command" error
-------------------------------------------------------------------------
1. Make sure Python is loaded and working:

.. code-block:: matlab

   >> py.print('test')

   test

2. Make sure Lentil is loaded and working:

.. code-block:: matlab

    >> mask = py.lentil.circle([int16(256),int16(256)],int16(128));

3. Verify there are no import errors in the Python code by importing Lentil and any
custom models in a Python interpreter:

.. code-block:: pycon

    >>> import lentil
    >>> import <<your-model>>

For more hints, see the MATLAB documentation on `Undefined variable "py" or function
"py.command" <https://www.mathworks.com/help/matlab/matlab_external/undefined-variable-py-or-function-py-command.html>`_

Resolving "Python Error: ImportError: Importing the numpy c-extensions failed." error
-------------------------------------------------------------------------------------
On Windows, if the system path is not correctly configured, Python will throw a lengthy
error message when trying to import Numpy:

.. code-block:: matlab

    >> py.importlib.import_module('numpy')
    Error using __init__><module> (line 54)
    Python Error: ImportError:

    IMPORTANT: PLEASE READ THIS FOR ADVICE ON HOW TO SOLVE THIS ISSUE!

    Importing the numpy c-extensions failed. - Try uninstalling and reinstalling numpy.

    ...

There appear to be several causes, but the error is most likely triggered because Python
is not able to locate the necessary Numpy DLL files. The most common culprit is not
electing to add Anaconda to the Windows PATH (which for some reason is the recommended
choice during installation). The issue is fixed by appending the system path. Note that
it is safest to do this from within MATLAB, in case a different version of Python is
on the system path and is being used by other applications.

.. code-block:: matlab

    setenv('path',['C:\Path\To\Anaconda3\Library\bin;', getenv('path')]);


Useful Links
============

* `MATLAB Examples <../examples>`_
* `Handling Data Returned from Python <https://www.mathworks.com/help/matlab/matlab_external/handling-data-returned-from-python.html>`_
* `Limitations to Python Support <https://www.mathworks.com/help/matlab/matlab_external/limitations-to-python-support.html>`_
* `Reloading Modified User-Defined Python Modules <https://www.mathworks.com/help/matlab/matlab_external/call-modified-python-module.html>`_
