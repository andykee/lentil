.. _user.advanced.extend:

****************
Extending Lentil
****************

Customizing Plane
=================
The Plane class or any of the classes derived from Plane can be subclassed to modify
any of the default behavior. Reasons to do this may include but are not limited to:

* Dynamically computing the :attr:`~lentil.Plane.opd` attribute
* Changing the Plane-Wavefront interaction by redefining the `Plane.multiply()` method
* Modifying the way a Plane is resampled or rescaled

Some general guidance for how to safely subclass Plane is provided below.

.. note::

    Lentil's |Plane| class and its subclasses all use Python's ``__init_subclass__()``
    method to ensure any required default values are set - even if a user-defined
    subclass does not explicitly call ``Plane``'s constructor ``__init__()`` method. For
    this reason, it is not strictly necessary to call ``super().__init__()`` when
    implementing a custom Plane subclass. It also won't hurt, as long as you're careful
    to either call ``super().__init__()`` before defining any static plane attributes or
    passing these attributes along to the ``super().__init__()`` call to ensure they are
    properly set.

Redefining the amplitude, OPD, or mask attributes
---------------------------------------------------
Plane :attr:`~lentil.Plane.amplitude`, :attr:`~lentil.Plane.opd`, and
:attr:`~lentil.Plane.mask` are all defined as properties, but Python allows you to
redefine them as class attributes without issue:

.. code-block:: python3

    import lentil

    class CustomPlane(le.Plane):
        def __init__(self):
            self.amplitude = lentil.circle((256,256), 128)
            self.opd = lentil.zernike(lentil.circlemask((256,256),128), 4)

If more dynamic behavior is required, the property can be redefined. For example, to
return a new random OPD each time the :attr:`~lentil.Plane.opd` attribute is
accessed:

.. code-block:: python3

    import numpy as np
    import lentil

    class CustomPlane(lentil.Plane):
        def __init__(self):
            self.mask = lentil.circlemask((256,256), 128)
            self.amplitude = lentil.circle((256,256), 128)

        @property
        def phaopdse(self):
            return lentil.zernike_compose(self.mask, np.random.random(10))

It is also straightforward to implement a custom :attr:`~lentil.Plane.opd` property to
provide a stateful OPD attribute:

.. code-block:: python3

    import numpy as np
    import lentil

    class CustomPlane(lentil.Plane):
        def __init__(self, x=np.zeros(10)):
            self.mask = lentil.circlemask((256,256), 128)
            self.amplitude = lentil.circle((256,256), 128)
            self.x = x

        @property
        def opd(self):
            return lentil.zernike_compose(self.mask, self.x)

.. note::

    Polychromatic or broadband diffraction propagations access the OPD, amplitude,
    and mask attributes for each propagatioon wavelength. Because these attributes
    remain fixed during a propagation, it is inefficient to repeatedly recompute
    them. To mitigate this, it can be very useful to provide a mechanism for freezing
    these dynamic attributes. There are many ways to do this. One approach is provided
    below:

    .. code-block:: python3

        import copy
        import numpy as np
        import lentil

        class CustomPlane(lentil.Plane):
            def __init__(self):
                self.mask = lentil.circlemask((256,256), 128)
                self.amplitude = lentil.circle((256,256), 128)

            @property
            def opd(self):
                return lentil.zernike_compose(self.mask, np.random.random(10))

            def freeze(self):
                # Return a copy of CustomPlane with the OPD attribute redefined
                # to be a static copy of the OPD when freeze() is called
                out = copy.deepcopy(self)
                out.opd = self.opd.copy()
                return out


Customizing Plane methods
-------------------------
Any of the |Plane| methods can be redefined in a subclass without restriction. Care
should be taken to ensure any redefined methods return data compatible with the
parent method's return type to preserve compatibility within Lentil.


Implementing the _Tilt interface
--------------------------------