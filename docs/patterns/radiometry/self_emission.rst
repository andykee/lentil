Including Self-Emission Flux in IR Models
=========================================

.. code:: python3

    class Model:

        ...

        def background_flux(self, qe=1):

            pixel_area = self.detector.pixelscale**2
            optics = pixel_area * self.pupil.omega * self.pupil.emission

            return optics * qe

        def image(...):

            ...


            # do the propagation
            img = lentil.propagate(...)

            # include the background
            background = self.background_flux(qe=self.detector.qe)
            background = background.integrate() * np.ones_like(img)
            background = background / oversample**2  # We have to account for the fact that our pixels are oversampled at this point
            img += background

            # simulate a frame
            frame = self.detector.frame(flux=img, ...)
