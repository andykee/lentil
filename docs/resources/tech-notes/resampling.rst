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
    sc = ax.scatter(cols.ravel(), rows.ravel(), c=psf_rebin.ravel()**0.5, s=20, marker='o', cmap='turbo')
    ax.set_aspect('equal')
    ax.invert_yaxis()

    ax.set_xticks(np.linspace(-0.53, psf_rebin.shape[1]-0.5, shape[1]+1))
    ax.set_yticks(np.linspace(-0.53, psf_rebin.shape[0]-0.5, shape[0]+1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(length=0)
    ax.grid(color='#333')

    ax.set_facecolor('#eee')