import matplotlib.pyplot as plt
import lentil

pupil = lentil.Pupil(amplitude=lentil.util.circle((256, 256), 128), diameter=1,
                     focal_length=10, pixelscale=1/256)

detector = lentil.Detector(pixelscale=5e-6, shape=(1024, 1024))
psf = lentil.propagate([pupil, detector], wave=650e-9, npix=(64, 64))
plt.imshow(psf, origin='lower')
plt.savefig('../../_static/img/psf_64.png', transparent=True, bbox_inches='tight', dpi=150)

tilt = lentil.Tilt(x=10e-6, y=-5e-6)
psf_tilt = lentil.propagate([tilt, pupil, detector], wave=650e-9, npix=(64, 64))
plt.imshow(psf_tilt, origin='lower')
plt.savefig('../../_static/img/psf_64_tilt.png', transparent=True, bbox_inches='tight', dpi=150)
