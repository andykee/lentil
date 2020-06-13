.. _cookbook-matlab:

Simple MATLAB Interface Class
=============================

Consider the following simple Lentil model which consists of an instrument, optical
system, and detector:

.. code-block:: python3

    import lentil

    class TinyInstrument(mo.Instrument):
        optical_system = TinyOpticalSystem()
        detector = TinyDetector()

    class TinyOpticalSystem(mo.SimpleOpticalSystem):
        diameter = 1
        focal_length = 8
        pixelscale = self.diameter/512
        detector_pixelscale = 5e-6

        def amplitude(self):
            return le.util.circle((512,512), 256)

    class TinyDetector(mo.FPA):
        pixelscale = 5e-6
        shape = (512,512)
        qe = le.radiometry.Spectrum(wave=[400, 600, 1000], value=[0.4, 0.8, 0.05],\
                         waveunits='nm', valueunits=None)
        gain = le.detector.Gain(gain=2.44, saturation_capacity=10000)

A corresponding MATLAB class that provides an interface to the instrument looks like
this:

.. code-block:: matlab

    classdef TinyInstrument < handle

        properties(Hidden):
            cls
        end

        methods
            function self = TinyInstrument
                self.cls = py.tiny.TinyInstrument()
            end
        end

    end

The MATLAB interface can be easily packaged up for delivery with the Python model by
including it in a to-level subdirectory:

.. code-block:: none
    :emphasize-lines: 12,13,14,15

    ~/dev/tiny-lentil
    ├── tiny_telescope/
    │   ├── __init__.py
    │   ├── detector.py
    │   ├── tiny.py
    │   ├── pupil.py
    │   ├── radiometry.py
    │   └── data/
    │       ├── detector_qe.csv
    │       ├── pupil_mask.npy
    │       └── sensitivities.npz
    ├── matlab/
    │   ├── mat2ndarray.m
    │   ├── ndarray2mat.m
    │   └── tiny.m
    ├── docs/
    ├── scripts/
    ├── tests/
    ├── .gitignore
    ├── README.md
    └── setup.py

Note the inclusion of ``mat2ndarray.m`` and ``ndarray2mat``.
