import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.circle(shape=(256,256), radius=120)
amp -= lentil.circle(shape=(256,256), radius=40)
for angle in (0, 90, 180, 270):
    amp *= lentil.spider((256,256), 2, angle)
opd = lentil.zernike_compose(amp, (0, 3.5e-6, -2e-6))
pupil = lentil.Pupil(amplitude=amp, opd=opd, pixelscale=1/240,
                     focal_length=10)
w0 = lentil.Wavefront(wavelength=500e-9)
w1 = w0 * pupil
w2 = lentil.propagate_dft(w1, shape=(64,64), pixelscale=5e-6, oversample=4)
fig, ax = plt.subplots()
ax.imshow(w2.intensity, norm='log', vmin=10e-5, cmap='inferno')
ax.axis('off')
fig.savefig('tilt.png', dpi=150, bbox_inches='tight')
