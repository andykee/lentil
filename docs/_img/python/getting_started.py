import numpy as np
import matplotlib.pyplot as plt
import monocle as mo

pm_mask = mo.util.circle(shape=(512, 512), radius=256)
sm_obsc = mo.util.circle(shape=(512, 512), radius=256/3)
sm_obsc_vert = np.zeros((512,512))
sm_obsc_vert[:, 252:258] = 1
sm_obsc_horiz = np.zeros((512, 512))
sm_obsc_horiz[252:258, :] = 1
mask = pm_mask - sm_obsc - sm_obsc_vert - sm_obsc_horiz
mask[np.where(mask < 0)] = 0  # Make the mask binary (0's and 1's)
plt.imshow(mask)
plt.savefig('../../_static/img/quickstart_pupil.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()

simple_pupil = mo.Pupil(focal_length=10, pixelscale=0.5/512, amplitude=mask)
simple_detector = mo.Detector(pixelscale=15e-6, shape=(512, 512))
planes = [simple_pupil, simple_detector]
psf = mo.propagate(planes, wave=550e-9, npix=(128,128))
plt.imshow(psf**0.1)
plt.savefig('../../_static/img/quickstart_psf_550_native.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()


psf = mo.propagate(planes, wave=np.arange(450e-9,650e-9,10e-9), npix=(128,128), oversample=5, rebin=False)

plt.imshow(psf**0.1)
plt.savefig('../../_static/img/quickstart_psf_broadband_5.png', transparent=True, bbox_inches='tight', dpi=150)


simple_pupil = mo.Pupil(focal_length=10,
                        pixelscale=0.5/512,
                        amplitude=mask,
                        phase=-1400e-9 * mo.zernike.zernike(mask, index=11))
plt.imshow(simple_pupil.phase)
plt.savefig('../../_static/img/quickstart_opd_spherical.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()

planes = [simple_pupil, simple_detector]
psf = mo.propagate(planes, wave=np.arange(475e-9,485e-9,1e-9), npix=(128,128))
plt.imshow(psf**0.1)
plt.savefig('../../_static/img/quickstart_psf_480_spherical.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()


qe = mo.radiometry.Spectrum(wave=[400, 600, 1000], value=[0.4, 0.8, 0.05], waveunit='nm', valueunit=None)
gain = mo.detector.Gain(gain=0.0016, saturation_capacity=10000)
simple_detector = mo.FPA(pixelscale=15e-6, shape=(512, 512), qe=qe, gain=gain)
irrad = np.tile(np.arange(512)/511, (512,1))  # gradient spanning [0, 1]
flux = irrad * 15000  # flux spanning [0, 15000] e-
plt.imshow(flux, cmap='gray')
plt.colorbar()
plt.savefig('../../_static/img/quickstart_detector_irradiance.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()


img = simple_detector.frame(flux, ts=1.5, wave=650, waveunit='nm')
plt.imshow(img, cmap='gray')
plt.colorbar()
plt.savefig('../../_static/img/quickstart_detector_frame.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()


detector_attrs = {
        'qe': qe,
        'gain': mo.detector.Gain(gain=2**12/15000, saturation_capacity=15000),
        'shot_noise': mo.detector.ShotNoise(),
        'read_noise': mo.detector.ReadNoise(50),
        'dark_signal': mo.detector.DarkCurrent(500)}

simple_detector = mo.FPA(pixelscale=15e-6, shape=(512, 512), **detector_attrs)
simple_pupil = mo.Pupil(focal_length=10,
                        pixelscale=0.5/512,
                        amplitude=mo.util.normalize_power(mask),
                        phase=-1400e-9 * mo.zernike.zernike(mask, index=11))

src = mo.radiometry.Blackbody.vegamag(wave=np.arange(350, 750),
                                      temp=9500,
                                      mag=8,
                                      band='V')

ir_filter = mo.radiometry.Spectrum(wave=[350, 395, 400, 700, 705, 750], value=[0, 0, 0.8, 0.8, 0, 0])
collecting_area = np.pi*(simple_pupil.diameter/2 - simple_pupil.diameter/6)**2
irradiance = src * ir_filter * collecting_area

fig, ax = plt.subplots(3,1, figsize=[5, 5])

ax[0].plot(src.wave, src.value)
ax[0].grid()
ax[0].set_ylim([40000,80000])
ax[0].set_title('Stellar flux')
ax[0].set_ylabel('ph/s/m^2/λ')

ax[1].plot(ir_filter.wave, ir_filter.value)
ax[1].grid()
ax[1].set_ylim([0,1])
ax[1].set_title('IR filter')
ax[1].set_ylabel('a.u.')

ax[2].plot(irradiance.wave, irradiance.value)
ax[2].grid()
ax[2].set_ylim([0,6000])
ax[2].set_title('Detector irradiance')
ax[2].set_ylabel('ph/s/λ')
ax[2].set_xlabel('Wavelength [nm]')

plt.tight_layout()
plt.savefig('../../_static/img/quickstart_source.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()


wave = np.arange(350,750,10)
binned_irradiance = irradiance.bin(wave)
planes = [simple_pupil, simple_detector]
psf = mo.propagate(planes, wave=wave*1e-9, weight=binned_irradiance,
                   npix=(128, 128), flatten=False)
img = simple_detector.frame(psf, ts=1e-4, wave=wave, collect_charge=True)
plt.imshow(img, cmap='gray')
plt.savefig('../../_static/img/quickstart_img.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()
