import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.circle(shape=(256,256), radius=120)
pupil = lentil.Pupil(amplitude=amp, opd=0, pixelscale=1/240,
                     focal_length=20)
w0 = lentil.Wavefront(wavelength=500e-9)
w1 = w0 * pupil
w2 = lentil.propagate_dft(w1, shape=(64,64), pixelscale=5e-6, oversample=4)
fig, ax = plt.subplots()
ax.imshow(w2.intensity, norm='log', vmin=10e-5)
ax.axis('off')
fig.savefig('simple_diffraction.png', dpi=150, bbox_inches='tight')
