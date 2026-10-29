Resampling oversampled images to native detector sampling
=========================================================


.. plot::
    :context: reset
    :scale: 50

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
    
    fig, ax = plt.subplots()
    ax.imshow(psf, cmap='turbo')
    ax.axis('off')

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(psf, cmap='turbo')
    ax.set_xticks(np.linspace(0, psf.shape[1], shape[1]+1))
    ax.set_yticks(np.linspace(0, psf.shape[0], shape[0]+1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(length=0)
    ax.grid(color='#999')

.. plot::
    :context: close-figs
    :scale: 50

    psf_rebin = lentil.rebin(psf, oversample)
    fig, ax = plt.subplots()
    rows, cols = np.indices(psf_rebin.shape)
    sc = ax.scatter(cols.ravel(), rows.ravel(), c=psf_rebin.ravel()**0.5, s=5, marker='o', cmap='turbo')
    ax.set_aspect('equal')
    ax.invert_yaxis()

    ax.set_xticks(np.linspace(-0.53, psf_rebin.shape[1]-0.5, shape[1]+1))
    ax.set_yticks(np.linspace(-0.53, psf_rebin.shape[0]-0.5, shape[0]+1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(length=0)
    ax.grid(color='#333')

    ax.set_facecolor('#eee')


.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(psf_rebin, cmap='turbo')
    ax.axis('off')


.. plot::
    :context: close-figs
    :scale: 50

    psf_discrete = lentil.sample(psf, 25)

    fig, ax = plt.subplots()
    rows, cols = np.indices(psf_discrete.shape)
    sc = ax.scatter(cols.ravel(), rows.ravel(), c=psf_discrete.ravel()**0.5, s=0.5, marker='o', cmap='turbo')
    ax.set_aspect('equal')
    
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(length=0)
    #ax.grid(color='#aaa')
    ax.set_xlim([-0.5,47.5])  # 59.9 when take = 20, 47.5 when take = 25
    ax.set_ylim([-0.5,47.5])

    ax.set_facecolor('#eee')

    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


.. plot::
    :context: close-figs
    :scale: 50

    psf_pixel = lentil.detector.pixel(psf, oversample=100)
    psf_pixel_native = lentil.sample(psf_pixel, 100)
    fig, ax = plt.subplots()
    ax.imshow(psf_pixel_native, cmap='turbo')
    ax.axis('off')


.. plot::
    :context: close-figs
    :scale: 50

    import itertools

    psf_pixel_nowrap = np.copy(psf_pixel)
    n = 285
    psf_pixel_nowrap[n:,:] = psf_pixel_nowrap[n-1,:]
    psf_pixel_nowrap = psf_pixel

    s = 15
    a = lentil.sample(psf_pixel_nowrap, oversample, offset=(-s,-s))
    b = lentil.sample(psf_pixel_nowrap, oversample, offset=(-s,0))
    c = lentil.sample(psf_pixel_nowrap, oversample, offset=(-s,s))
    d = lentil.sample(psf_pixel_nowrap, oversample, offset=(0,-s))
    e = lentil.sample(psf_pixel_nowrap, oversample, offset=(0,0))
    f = lentil.sample(psf_pixel_nowrap, oversample, offset=(0,s))
    g = lentil.sample(psf_pixel_nowrap, oversample, offset=(s,-s))
    h = lentil.sample(psf_pixel_nowrap, oversample, offset=(s,0))
    i = lentil.sample(psf_pixel_nowrap, oversample, offset=(s,s))

    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(5,5))
    ax[0,0].imshow(a, cmap='turbo')
    ax[0,1].imshow(b, cmap='turbo')
    ax[0,2].imshow(c, cmap='turbo')
    ax[1,0].imshow(d, cmap='turbo')
    ax[1,1].imshow(e, cmap='turbo')
    ax[1,2].imshow(f, cmap='turbo')
    ax[2,0].imshow(g, cmap='turbo')
    ax[2,1].imshow(h, cmap='turbo')
    ax[2,2].imshow(i, cmap='turbo')

    for a in itertools.product(range(3), range(3)):
        ax[a].axis('off')