.. _examples.bandpass_resampling:

****************************
Bandpass and flux resampling
****************************

Resampling a bandpass
---------------------


Bandpass represented by arrays:

.. code:: python

    wave = np.array([...])
    bandpass = np.array([...])

    wave_sampling = 5e-9 # wavelength sampling interval
    
    
    start, stop = wave[0], wave[-1]
    num = int(round(sto-start) / (wave_sampling*1e9))

    resampled_wave = np.linspace(start, stop, num)
    resampled_bandpass = np.interp(resampled_wave, wave, bandpass)




Bandpass represented by a :class:`~lentil.radiometry.Spectrum` object:

.. code:: python

    bandpass = lentil.Spectrum(...)

    wave_sampling = 5e-9 # wavelength sampling interval
    trim_tol = 1e-2 # relative tolerance used to clip off nearly zero ends

    bandpass.trim(trim_tol)
    start, stop = flux.wave[0], flux.wave[-1]
    num = int(round(sto-start) / wave_sampling*1e9)

    resampled_wave = np.linspace(start, stop, num)
    resampled_bandpass = bandpass.sample(wave, waveunit=bandpass.waveunit)



Resampling and rebinning flux
-----------------------------
In this example, the flux is binned (instead of simply sampled) to preserve
its integrated power.

.. code:: python

    flux = lentil.Spectrum(...)

    wave_sampling = 5e-9 # wavelength sampling interval
    trim_tol = 1e-2 # relative tolerance used to clip off nearly zero ends

    flux.trim(ftrim_tol)
    start, stop = flux.wave[0], flux.wave[-1]
    num = int(round(sto-start) / (wave_sampling*1e9))

    wave = np.linspace(start, stop, num)
    binned_flux = flux.bin(wave, waveunit=flux.waveunit)
