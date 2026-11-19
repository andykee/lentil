Resampling oversampled images to native detector sampling
=========================================================

Diffraction calculations are usually performed on a grid that is finer than the
detector's pixel pitch. Oversampling ensures the intensity is not aliased (the
optical cutoff frequency of a typical imaging system lies above the detector's
Nyquist frequency) and makes it possible to represent effects that occur on
sub-pixel scales. This note describes two equivalent approaches for resampling
oversampled data to a detector's native sampling.

Pixel integration is convolution followed by sampling
-----------------------------------------------------

Before we describe the two resampling approaches it is worth exploring the
spatial aspect of how a detector digitizes a light field. An ideal square pixel
with pitch :math:`p` and a 100% fill factor collects all of the power falling
on its area (integration). For a detector with :math:`N_r` rows and
:math:`N_c` columns, the signal in the pixel in row :math:`r` and column
:math:`c` (:math:`r = 0 \ldots N_r - 1`, :math:`c = 0 \ldots N_c - 1`) is

.. math::
    :label:

    P[r,c] = \iint_{-\infty}^{\infty} I(u,v) \,
             \Pi\left(\frac{u-u_c}{p}\right) \Pi\left(\frac{v-v_r}{p}\right)
             \, du \, dv \, ,

where :math:`I(u,v)` is the irradiance in the image plane, :math:`\Pi` is the
rectangle function, and :math:`(u_c, v_r)` is the center of the pixel.
Equation 1 is a two-dimensional convolution of the image plane intensity
with the pixel aperture function. Defining the pixel aperture function as
:math:`\Pi_{px} = \Pi(u/p)\Pi(v/p)`, equation 1 becomes

.. math::
    :label:

    P[r,c] = \left(I * \Pi_{px}\right)(u_c, v_r) \, .

Integrating the image over each pixel is exactly equivalent to convolving the
image with the pixel aperture and then point-sampling the result at the
pixel centers. The rest of this note will show how to implement both of these
approaches using built-in lentil functions. Let :math:`f` be the oversampled
PSF computed with an oversampling factor of :math:`M = 100` (Figure 1(a)) on
a :math:`N_r \times N_c = 12 \times 12` detector grid with pixels spaced
:math:`p/M` (Figure 1 (b)).

.. plot::
    :context: reset
    :width: 70%
    :caption: **Figure 1**: Oversampled PSF

    amp = lentil.circle((256,256), 120)
    opd = lentil.zernike_compose(amp, [0,0,0,10e-9,-45e-9,12e-9,-75e-9,8e-9,0,-16e-9])
    p = lentil.Pupil(amplitude=amp, opd=opd, mask=amp, pixelscale=1/240,
                     focal_length=9, diameter=1)
    shape = (12,12)
    oversample = 100
    psf = 0
    for wl in [450e-9, 500e-9, 550e-9]:
        w = lentil.Wavefront(wavelength=wl)
        w = w * p
        w = lentil.propagate_dft(w, shape=shape, pixelscale=5e-6, oversample=oversample)
        psf += w.intensity

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
    ax[0].imshow(psf, cmap='turbo')
    ax[0].set_xlabel('(a)')
    ax[0].set_xticklabels([])
    ax[0].set_yticklabels([])
    ax[0].tick_params(length=0)

    ax[1].imshow(psf, cmap='turbo')
    ax[1].set_xticks(np.linspace(0, psf.shape[1], shape[1]+1))
    ax[1].set_yticks(np.linspace(0, psf.shape[0], shape[0]+1))
    ax[1].set_xticklabels([])
    ax[1].set_yticklabels([])
    ax[1].tick_params(length=0)
    ax[1].grid(color='#999')
    ax[1].set_xlabel('(b)')



Method 1: Rebinning
-------------------

The direct (rebinning) approach evaluates the integral in Equation 1 by summing
the :math:`M \times M` block of oversampled samples that fall inside each
pixel:

.. math::
    :label:

    P[r,c] = \sum_{i=0}^{M-1} \sum_{j=0}^{M-1} f[rM + i, \, cM + j] \, ,

where :math:`i` and :math:`j` index the sub-pixel samples within each block.
This functionality is provided by lentil's :func:`~lentil.rebin` function:

.. code:: pycon

   >>> psf_rebin = lentil.rebin(psf, factor=100)

.. plot::
    :context: close-figs
    :scale: 50
    :caption: **Figure 2**: Rebinned PSF at native detector sampling.

    psf_rebin = lentil.rebin(psf, oversample)
    fig, ax = plt.subplots()
    ax.imshow(psf_rebin, cmap='turbo')
    ax.axis('off')


Method 2: Pixel MTF convolution and sampling
--------------------------------------------

The second approach implements Equation 2 by first convolving the oversampled
image with the pixel aperture and then extracting one sample per pixel.

**1. Convolution:** For a square pixel with 100% fill factor, its aperture
function is a :math:`\mathrm{rect}()` function in the spatial domain or a
:math:`\mathrm{sinc}()` function in the frequency domain. By the convolution
theorem, we can evaluate Equation 2 by multiplying the image spectrum with
the frequency domain representation of the pixel aperture function:

.. math::
    :label:

    g = \mathcal{F}^{-1}\left\{ \mathcal{F}\{f\} \cdot
        \mathrm{sinc}(M\nu_x) \, \mathrm{sinc}(M\nu_y) \right\} \, ,

where :math:`\nu_x, \nu_y` are frequencies in cycles per oversampled image
sample.

**2. Sample:** Select one sample from each oversampled :math:`M \times M`
block:

.. math::
    :label:

    s[r,c] = g[rM + q_y, \, cM + q_x] \, , \qquad q = \lfloor M/2 \rfloor + \mathrm{offset} \, ,

where ``offset`` selects the sub-pixel phasing (details below). Note that
because :math:`g` holds pixel averages of :math:`f`, the result must be
multiplied by :math:`M^2` to recover the power in each pixel:

.. math::
    :label:

    P[r,c] = M^2 \, s[r,c] \, .

These operations are implemented in lentil via the :func:`lentil.detector.pixel` and
:func:`~lentil.sample` functions:

.. code:: pycon

    >>> psf_pixel = lentil.detector.pixel(psf, oversample=100)
    >>> psf_pixel_sample = lentil.sample(psf_pixel, oversample=oversample)
    >>> psf_pixel_sample *= 100**2

.. plot::
    :context: close-figs
    :scale: 50
    :caption: **Figure 3**: PSF convolved with the pixel MTF and then sampled at
              native detector sampling.

    psf_pixel = lentil.detector.pixel(psf, oversample=oversample)
    psf_pixel_native = lentil.sample(psf_pixel, oversample) * oversample**2
    fig, ax = plt.subplots()
    ax.imshow(psf_pixel_native, cmap='turbo')
    ax.axis('off')


Pixel phasing
-------------

The ``offset`` argument to :func:`lentil.sample` chooses which sample within
each block is kept. This is equivalent to shifting the pixel grid relative to
the image in steps of :math:`1/M` pixel. An offset of :math:`+o` samples moves
the image by :math:`-o/M` pixels in the resampled result. Figure 7 shows the
natively sampled PSF for offsets of :math:`\pm 15` samples (:math:`\pm 0.15`
pixel) in each axis. All nine images are computed from the same convolved
array :math:`g`, and each contains the same total power. This makes method 2 a
convenient way to generate PSFs at many sub-pixel positions, for example to
model pointing jitter or dithered observations, from a single diffraction
calculation.

.. plot::
    :context: close-figs
    :scale: 50
    :caption: **Figure 7**: Natively sampled PSF for nine pixel phasings. The
              offset in :math:`y` varies down the rows and the offset in
              :math:`x` varies across the columns, in steps of 0.15 pixel.

    import itertools

    s = 15
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(5,5))
    for (i, dy), (j, dx) in itertools.product(enumerate((-s, 0, s)), repeat=2):
        img = lentil.sample(psf_pixel, oversample, offset=(dx, dy)) * oversample**2
        ax[i,j].imshow(img, cmap='turbo')
        ax[i,j].axis('off')

Rebinning can produce exactly the same phasings by shifting the block
boundaries, i.e. by cropping the oversampled array :math:`o` samples later
before calling :func:`lentil.rebin`.


Edge effects
------------

The Fourier transform in Equation 4 performs a *circular* convolution so values
near an edge of the array wrap around to the opposite edge. With ``offset=0``
and odd :math:`M` the pixel-sized window stays within its own block, so only
the small ringing of the band-limited kernel wraps. For any other phasing, the
windows of the edge pixels extend past the array and pick up samples from the
far side. To avoid this effect, compute the oversampled PSF over a field at
least one native pixel larger than needed on each side, resample, and then crop
the edge pixels.

