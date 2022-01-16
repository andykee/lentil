import numpy as np
import lentil as le
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = (6.5, 6.5)

diameter = 1
focal_length = 5
du = 5e-6
wl = 600e-9
npix_pupil = 256
npix_image = 32
oversample = np.array([1,2,3,4])
dx = diameter/npix_pupil

amp = le.util.circle((npix_pupil,npix_pupil), npix_pupil//2)
f_number = focal_length/diameter

#Q = (wl*f_number)/du
#Qd = Q*oversample
#q = 1/(alpha * npix)
#print(f'System Q: {Q}       Discrete Q: {Qd}       q: {q}')

psf = []
for os in oversample:
    npix = npix_image * os
    alpha = (dx*du)/(wl*focal_length*os)
    F = le.fourier.dft2(amp, alpha, npix)
    p = le.util.rebin(np.abs(F)**2, os)
    psf.append(p)

exp = 0.2

plt.subplot(1,4,1)
plt.imshow(psf[0]**exp)
plt.title('$Q = 0.6$')
plt.axis('off')

plt.subplot(1,4,2)
plt.imshow(psf[1]**exp)
plt.title('$Q = 1.2$')
plt.axis('off')

plt.subplot(1,4,3)
plt.imshow(psf[2]**exp)
plt.title('$Q = 1.8$')
plt.axis('off')

plt.subplot(1,4,4)
plt.imshow(psf[3]**exp)
plt.title('$Q = 2.4$')
plt.axis('off')
