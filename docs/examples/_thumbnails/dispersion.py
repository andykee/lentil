import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.circle(shape=(256,256), radius=120)

pupil = lentil.Pupil(amplitude=amp, opd=0, pixelscale=1/240,
                     focal_length=10)

trace = [1, 0] 
dispersion = [3e-4, 550e-9]
dispersive_element = lentil.DispersiveTilt(trace=trace, dispersion=dispersion)


shape = (128,128)
oversample = 3
out = np.zeros((shape[0]*oversample, shape[1]*oversample))

for wave in np.linspace(475, 625, 150):
    w0 = lentil.Wavefront(wavelength=wave*1e-9)
    w1 = w0 * pupil
    w2 = w1 * dispersive_element
    w3 = lentil.propagate_dft(w2, shape=shape, prop_shape=(32, 32), 
                              pixelscale=5.5e-6, oversample=oversample)

    out = w3.insert(out)

fig, ax = plt.subplots()
ax.imshow(out, cmap='inferno')
ax.axis('off')
fig.savefig('dispersion.png', dpi=150, bbox_inches='tight')