import matplotlib.pyplot as plt
import lentil

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4.5, 4.5)

mask = lentil.circle((256, 256), 120)
opd = lentil.zernike_compose(mask, [0, 0, 0, -300e-9, 50e-9, -100e-9, 50e-9])
pupil = lentil.Pupil(amplitude=mask, opd=opd, focal_length=10, pixelscale=1 / 240)
w = lentil.Wavefront(650e-9)
w *= pupil
w = lentil.propagate_dft(w, pixelscale=5e-6, shape=64, oversample=5)
psf = w.intensity

psf_smear = lentil.smear(psf, distance=5e-5,
                         pixelscale=5e-6, oversample=5)
plt.subplot(121)
plt.imshow(psf, cmap='inferno')
plt.title('Input image')
plt.subplot(122)
plt.imshow(psf_smear, cmap='inferno')
plt.title('Smeared image')
