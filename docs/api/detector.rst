.. _api-detector:

lentil.detector
===============

Gain
----
.. autosummary::

    lentil.detector.Gain
    lentil.detector.PolynomialGain

Dark Signals
------------
.. autosummary::

    lentil.detector.DarkSignal
    lentil.detector.DarkCurrent
    lentil.detector.Rule07DarkCurrent

Random Noise
------------
.. autosummary::

    lentil.detector.RandomNoise
    lentil.detector.ShotNoise
    lentil.detector.GaussianShotNoise
    lentil.detector.ReadNoise
    lentil.detector.ChargeDiffusion

Fixed Pattern Noise
-------------------
.. autosummary::

    lentil.detector.FPN
    lentil.detector.NormalFPN
    lentil.detector.LognormalFPN
    lentil.detector.PRNU
    lentil.detector.DarkCurrentFPN
    lentil.detector.PixelOffsetFPN
    lentil.detector.ColumnOffsetFPN

Defects
-------
.. autosummary::

    lentil.detector.Defect
    lentil.detector.PixelMask
    lentil.detector.CosmicRay


.. Gain

.. autoclass:: lentil.detector.Gain
    :members:
    :special-members: __call__

.. autoclass:: lentil.detector.PolynomialGain
    :members:
    :inherited-members:
    :special-members: __call__

.. Dark Signals

.. autoclass:: lentil.detector.DarkSignal
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.DarkCurrent
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.Rule07DarkCurrent
    :members:
    :inherited-members:
    :special-members: __call__

.. Random noise

.. autoclass:: lentil.detector.RandomNoise
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.ShotNoise
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.GaussianShotNoise
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.ReadNoise
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.ChargeDiffusion
    :members:
    :inherited-members:
    :special-members: __call__

.. Fixed Pattern Noise

.. autoclass:: lentil.detector.FPN
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.NormalFPN
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.LognormalFPN
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.PRNU
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.DarkCurrentFPN
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.PixelOffsetFPN
    :members:
    :inherited-members:
    :special-members: __call__

.. autoclass:: lentil.detector.ColumnOffsetFPN
    :members:
    :inherited-members:
    :special-members: __call__

.. Defects

.. autoclass:: lentil.detector.Defect
    :members:
    :special-members: __call__

.. autoclass:: lentil.detector.PixelMask
    :members:
    :special-members: __call__

.. autoclass:: lentil.detector.CosmicRay
    :members:
    :inherited-members:
    :special-members: __call__
