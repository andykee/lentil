import numpy as np
import lentil as le
import matplotlib.pyplot as plt

diameter = 1
focal_length = 10
du = 5e-6
wl = 650e-9
npix_pupil = np.array([16, 32, 64, 128]) # [q < 0.5, q < 1, 1 < q < 2, q > 2]
npix_image = 64
oversample = 5
npix = npix_image * oversample

f_number = focal_length/diameter


Q = (wl*f_number)/du
Qd = Q*oversample
dx = diameter/npix_pupil
alpha = (dx*du)/(wl*focal_length*oversample)
q = 1/(alpha * npix)

psf = []
for n in npix_pupil:
    dx = diameter/n
    alpha = (dx*du)/(wl*focal_length*oversample)
    amp = le.util.circle((n,n), n//2)
    F = le.fourier.dft2(amp, alpha, npix)
    psf.append(np.abs(F)**2)

exp = 0.15

plt.subplot(1,4,1)
plt.imshow(psf[0]**exp)
plt.title('$q = 0.325$', fontsize=9)
plt.axis('off')

plt.subplot(1,4,2)
plt.imshow(psf[1]**exp)
plt.title('$q = 0.65$', fontsize=9)
plt.axis('off')

plt.subplot(1,4,3)
plt.imshow(psf[2]**exp)
plt.title('$q = 1.3$', fontsize=9)
plt.axis('off')

plt.subplot(1,4,4)
plt.imshow(psf[3]**exp)
plt.title('$q = 2.6$', fontsize=9)
plt.axis('off')

plt.savefig('q_sweep.png', transparent=True, bbox_inches='tight', dpi=300)