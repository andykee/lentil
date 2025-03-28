.. _examples.matlab_interface:

MATLAB interface class
======================

Consider the following simple Lentil model which consists of an instrument, optical
system, and detector:

.. code-block:: python3

    import lentil

    class TinyInstrument:
        def __init__(self):
            self.pupil = lentil.Pupil(amplitude=lentil.util.circle((256,256), 128),
                                      focal_length=10, pixelscale=1/256)
            self.detector = lentil.Image(pixelscale=5e-6)

        @property
        def planes(self):
            return [self.pupil, self.detector]

        def propagate(self, wave, weight, npix=1024, oversample=3, rebin=True,
                      tilt='phase', npix_chip=None):
            return lentil.propagate(self.planes, wave, weight, npix, npix_chip,
                                    oversampe, rebin, tilt, flatten=True)

A corresponding MATLAB class that provides an interface to the instrument looks like
this:

.. code-block:: matlab

    classdef TinyInstrument < handle

        properties(Hidden)
            cls
        end

        methods
            function self = TinyInstrument
                self.cls = py.tiny.TinyInstrument()
            end

            function img = propagate(self, wave, weight, npix, oversample, rebin, tilt, npix_chip)
                if nargin < 8; npix_chip = py.None; end
                if nargin < 7; tilt = 'phase'; end
                if nargin < 6; rebin = true; end
                if nargin < 5; oversample = 3; end
                if nargin < 4; npix = 1024; end

                img = self.cls.propagate(wave, weight, int(npix), int(oversample),... 
                                         rebin, tilt, int(npix_chip));
            end

        end
    end

More recent versions of MATLAB (>= 2019b) support optional positional arguments. The
``propagate`` function can be rewritten to take advantage of this new functionality:

.. code-block:: matlab

    function img = propagate(self, wave, weight, npix, oversample, rebin, tilt, npix_chip)
        arguments
            self
            wave
            weight
            npix = 1024
            oversample = 3
            rebin = true
            tilt = 'phase'
            npix_chip = py.None
        end

        img = self.cls.propagate(wave, weight), int(npix), int(oversample),... 
                                 rebin, tilt, int(npix_chip));
    end


The MATLAB interface can be easily packaged up for delivery with the Python model by
including it in a to-level subdirectory:

.. code-block:: none
    :emphasize-lines: 12,13,14,15

    ~/dev/tiny-lentil
    ├── tiny_telescope/
    │   ├── __init__.py
    │   ├── planes.py
    │   ├── radiometry.py
    │   ├── telescope.py
    │   └── data/
    │       ├── dwdx.npy
    │       ├── detector_qe.csv
    │       └── pupil_mask.npy
    ├── matlab/
    │   └── tiny.m
    ├── docs/
    ├── scripts/
    ├── tests/
    ├── .gitignore
    ├── README.md
    └── setup.py


