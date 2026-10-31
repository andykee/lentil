import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import numpy as np
import lentil

# f/10, 650 nm, 5 um pixels -> Q = 1.3. Same 2.56 mm field of view, two
# different pupil samplings
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4, 2.3))
for a, n in zip(ax, (32, 128)):
    amp = lentil.circle((n, n), n//2)
    pupil = lentil.Pupil(amplitude=amp, pixelscale=1/n, focal_length=10)
    w = lentil.Wavefront(650e-9) * pupil
    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=512)
    a.imshow(w.intensity, cmap='turbo', norm=AsinhNorm(1e-4), origin='lower',
             extent=[-1.28, 1.28, -1.28, 1.28])
    period = 650e-9 * 10 * n * 1e3   # lambda z / dx [mm]
    a.set_title(f'{n} px across pupil\nperiod = {period:.2f} mm', fontsize=7)
    a.set_xlabel('[mm]')
    for k in (-1, 1):
        a.axvline(k*period/2, color='w', lw=0.5, ls='--')
        a.axhline(k*period/2, color='w', lw=0.5, ls='--')
    a.set_xlim(-1.28, 1.28)
    a.set_ylim(-1.28, 1.28)
fig.tight_layout()
