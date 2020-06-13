Implementing a Filter Wheel
===========================

Implementing a filter wheel requires managing state (we need some way to keep track of
which filter has been selected and is "active"). In most simple cases, changing filters
only impacts the radiometric properties of a model. More complicated cases like when
optical aberrations are filter dependent or when the filter wheel includes a grism
require changing both the radiometric and optical models.

Filter wheel function for managing radiometric model
----------------------------------------------------
When selecting different filters only changes the radiometric properties of a model,
a filter wheel can be represented by a simple ``bandpass`` function. The function
accepts a filter name as its only argument, and returns the model's optical bandpass:

.. code-block:: python3

    import numpy as np
    import lentil

    telescope = lentil.radiometry.Spectrum(wave=np.arange(400,751),
                                           value=0.8*np.ones(351),
                                           waveunit='nm',
                                           valueunit=None)

    open = lentil.radiometry.Spectrum(wave=np.arange(400,751),
                                      value=np.ones(351),
                                      waveunit='nm',
                                      valueunit=None)

    f420n = lentil.radiometry.Spectrum(wave=np.arange(410,431),
                                       value=np.ones(21),
                                       waveunit='nm',
                                       valueunit=None)

    def bandpass(filter_name):
        if filter.upper() == 'OPEN':
            return telescope * open
        elif filter.upper() == 'F420N':
            return telescope * f420n
        else:
            raise ValueError('Unknown filter')


Filter wheel function and custom code for managing radiometric and optical models
---------------------------------------------------------------------------------
The following approach is useful when a filter wheel selection also has an impact on the
optical model. In addition to updating the radiometric model, the reference to a filter
plane is adjusted in the model's ``planes`` list used for propagation. In this case,
we'll add a grism. First, we'll include the grism in the radiometric model:

.. code-block:: python3

    import numpy as np
    import lentil

    telescope = lentil.radiometry.Spectrum(wave=np.arange(400,751),
                                           value=0.8*np.ones(351),
                                           waveunit='nm',
                                           valueunit=None)

    open = lentil.radiometry.Spectrum(wave=np.arange(400,751),
                                      value=np.ones(351),
                                      waveunit='nm',
                                      valueunit=None)

    f420n = lentil.radiometry.Spectrum(wave=np.arange(410,431),
                                       value=np.ones(21),
                                       waveunit='nm',
                                       valueunit=None)

    grism = lentil.radiometry.Spectrum(wave=np.arange(450,651),
                                       value=np.ones(201),
                                       waveunit='nm',
                                       valueunit=None)


    def bandpass(filter_name):
        if filter.upper() == 'OPEN':
            return telescope * open
        elif filter.upper() == 'F420N':
            return telescope * f420n
        elif filter.upper() == 'GRISM':
            return telescope * grism
        else:
            raise ValueError('Unknown filter')


Next, we'll import the radiometric model and update our optical model to slot in the
``Grism`` plane when it is selected. Notice that because the other filters don't have
any optical effect, they are represented by clean ``Plane()`` objects.

.. code-block:: python3

    from examplemodel.planes import Pupil, Grism
    from examplemodel.detector import Detector
    from examplemodel import radiometry

    class ExampleModel:
        def __init__(self):
            self._pupil = Pupil()
            self._detector = Detector()
            self._filters = {'OPEN': mo.Plane(),
                             'F420N': mo.Plane(),
                             'GRISM': Grism()}

            self.filter = 'OPEN'

        @property
        def _planes(self):
            return [self._pupil, self._filter, self._detector]

        @property
        def filter(self):
            return self._filter_name

        @filter.setter
        def filter(self, name):
            try:
                self._filter = self._filters[name.upper()]
                self._filter_name = name.upper()
            except KeyError:
                raise ValueError('Unknown filter', name)

        @property
        def bandpass(self):
            return radiometry.bandpass(self.filter)
