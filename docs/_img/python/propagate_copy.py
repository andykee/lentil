import matplotlib.pyplot as plt
import lentil

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4.5, 4.5)

pupil = lentil.Pupil(amplitude=lentil.circle((256, 256), 120),
                     pixelscale=1 / 240, focal_length=10)

w1 = lentil.Wavefront(650e-9)
w2 = w1 * pupil
w3 = w2.propagate_image(pixelscale=5e-6, npix=64, oversample=5, inplace=False)

plt.subplot(121)
plt.imshow(w2.intensity, origin='lower')
plt.title('w2 intensity')

plt.subplot(122)
plt.imshow(w3.intensity ** 0.1, origin='lower')
plt.title('w3 intensity')
