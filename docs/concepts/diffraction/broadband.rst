Broadband propagation
=====================
A |Wavefront| is monochromatic, so to model a broadband source, propagate one
wavefront per wavelength and add the intensities:

.. plot::
    :context: reset
    :include-source:
    :nofigs:

    pupil = lentil.Pupil(amplitude=lentil.circle((256, 256), 120),
                         pixelscale=1/240, focal_length=10)

    wavelengths = np.arange(450, 650, 10)*1e-9
    img = np.zeros((320, 320))

    for wl in wavelengths:
        w = lentil.Wavefront(wl) * pupil
        w = lentil.propagate_dft(w, pixelscale=5e-6, shape=64, oversample=5)
        w.insert(img)

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(img, cmap=psf_cmap, norm=AsinhNorm(),
              extent=[-0.16, 0.16, -0.16, 0.16])
    ax.set_xlabel('[mm]')
    ax.set_ylabel('[mm]')

:func:`~lentil.Wavefront.insert` adds the wavefront's intensity directly into
``img``, without building a full intensity array for each wavelength. Pass
``weight`` to apply a spectral weight. ``img`` must have shape
``shape * oversample``.

Sampling across a band
----------------------
:math:`Q` and the wraparound period both shrink with wavelength, so check
sampling at the shortest wavelength in the band. See :doc:`sampling` for more.

Performance
-----------
* Use the smallest ``prop_shape`` that holds the PSF. See
  :doc:`output-plane`.
* Freeze planes whose OPD is expensive to compute so it is computed once
  instead of once per wavelength.
* For large output planes, :func:`~lentil.propagate_fft` with a
  preallocated ``scratch`` array avoids reallocating the padded array for
  every wavelength. Size it with :func:`~lentil.scratch_shape`.

See :doc:`/guides/performance` for more.