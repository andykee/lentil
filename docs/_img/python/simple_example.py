import matplotlib.pyplot as plt
import lentil

mask = lentil.util.circle((256, 256), 128) - lentil.util.circle((256, 256), 128/3)
opd = lentil.zernike.zernike_compose(mask, coeffs=[0, 0, 0, 300e-9, 50e-9, -100e-9, 50e-9])

pupil = lentil.Pupil(amplitude=mask, phase=opd, focal_length=10,
                     pixelscale=1/256)
detector = lentil.Image(pixelscale=5e-6)

psf5 = lentil.propagate([pupil, detector], wave=650e-9, npix=32, oversample=10, rebin=False)

psf = lentil.detector.pixelate(psf5, oversample=10)

plt.imshow(mask)
plt.savefig('../../_static/img/getting_started_amp.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()

plt.imshow(opd)
plt.savefig('../../_static/img/getting_started_opd.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()

plt.imshow(psf5)
plt.savefig('../../_static/img/getting_started_psf_oversample.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()

plt.imshow(psf)
plt.savefig('../../_static/img/getting_started_psf_detector.png', transparent=True, bbox_inches='tight', dpi=150)
plt.close()
