import matplotlib.pyplot as plt
import lentil

mask = lentil.util.circle((256, 256), 128)
opd = lentil.zernike.zernike(mask, 8) * 500e-9  # 500 nm of coma
pupil = lentil.Pupil(amplitude=mask, phase=opd, diameter=1,
                     focal_length=10, pixelscale=1/256)
rotation = lentil.Rotate(angle=30, unit='degrees')
flip = lentil.Flip(1)
detector = lentil.Detector(pixelscale=5e-6, shape=(1024, 1024))

psf = lentil.propagate([pupil, detector], wave=650e-9, npix=(128, 128))
plt.imshow(psf, origin='lower')
plt.savefig('../../_static/img/psf_coma.png', transparent=True, bbox_inches='tight', dpi=150)

psf = lentil.propagate([pupil, rotation, detector], wave=650e-9, npix=(128, 128))
plt.imshow(psf, origin='lower')
plt.savefig('../../_static/img/psf_coma_rotate.png', transparent=True, bbox_inches='tight', dpi=150)

psf = lentil.propagate([pupil, flip, detector], wave=650e-9, npix=(128, 128))
plt.imshow(psf, origin='lower')
plt.savefig('../../_static/img/psf_coma_flip.png', transparent=True, bbox_inches='tight', dpi=150)
