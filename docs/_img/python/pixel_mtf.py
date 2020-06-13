import numpy as np
import matplotlib.pyplot as plt
import monocle as mo

# Define some system parameters
diameter = 1
focal_length = 5
f_number = focal_length/diameter
du = 5e-6
wave = 550e-9

cutoff_freq = 1/(f_number*wave)  # cycles/m
cutoff_freq_mm = cutoff_freq/1e3  # cycles/mm

# Define some simulation parameters
npix = 256
oversample = 4
psf_sampling = du/oversample
psf_sampling_mm = psf_sampling * 1e3

# Set up the spatial frequencies we'll evaluate the MTF at
# Note that we would like to at least cover the cutoff frequency here
mtf_axis_mm = 1/psf_sampling_mm * np.arange(0,npix//2) / npix  # cycles/mm
assert mtf_axis_mm[-1] >= cutoff_freq_mm

mtf_axis_px = mtf_axis_mm * du/1e-3  # cycles/px

# Compute the analytical optical MTF
s = (mtf_axis_px*wave*f_number)/du
s_valid = s[s<1]
mtf_valid = 2/np.pi*(np.arccos(s_valid)-s_valid*np.sqrt(1-s_valid**2))
mtf_optics = np.zeros_like(s)
mtf_optics[0:s_valid.shape[0]] = mtf_valid

# Compute the analytical pixel MTF
f_px = mtf_axis_px
f_px[0] = 1e-15  # We have to cheat a little to avoid a divide
mtf_px = np.abs(np.sin(np.pi*f_px)/(np.pi*f_px))

# The system MTF is just the product of the optical and pixel MTFs
mtf_sys = mtf_optics * mtf_px

# Let's take a look
plt.style.use('ggplot')
#plt.plot(mtf_axis_px, mtf_optics, label='optics')
#plt.plot(mtf_axis_px, mtf_px, label='pixel')
#plt.plot(mtf_axis_px, mtf_sys, label='system')
#plt.xlabel('cycles/px')
#plt.ylabel('MTF')
#plt.legend()

npix_pup_rad = npix/2
npix_pup_diam = npix_pup_rad * 2
dx = diameter/npix_pup_diam
amp = mo.util.circle((npix,npix),npix_pup_rad)
alpha = (dx*du)/(wave*focal_length*oversample)

# Compute the optical MTF from a Monocle-generated PSF
psf = np.abs(mo.fourier.dft2(amp, alpha, npix=npix)**2)
psf = psf/np.max(psf)
mtf_optics_mo = np.abs(np.fft.fft2(psf))
mtf_optics_mo = mtf_optics_mo/np.max(mtf_optics_mo)
mtf_optics_mo = mtf_optics_mo[0,0:mtf_optics_mo.shape[0]//2]

# Now apply Monocle's pixellation method and compute the system MTF
pixel_mtf = mo.convolvable.Pixel()
psf_px = pixel_mtf(psf, oversample=oversample)
mtf_sys_mo = np.abs(np.fft.fft2(psf_px))
mtf_sys_mo = mtf_sys_mo/np.max(mtf_sys_mo)
mtf_sys_mo = mtf_sys_mo[0,0:mtf_sys_mo.shape[0]//2]

# Finally, we'll grab the Pixel kernel to make sure it matches the
# analytic pixel MTF
mtf_px_mo = np.abs(pixel_mtf.kernel(mtf_axis_px, mtf_axis_px, 1))
mtf_px_mo = mtf_px_mo[0,:]

#plt.plot(mtf_axis_px, mtf_optics_mo, label='optics')
#plt.plot(mtf_axis_px, mtf_px_mo, label='pixel')
#plt.plot(mtf_axis_px, mtf_sys_mo, label='system')
#plt.title('Model-based MTF')
#plt.ylabel('MTF')
#plt.xlabel('cycles/px')
#plt.legend()

# Finally we'll bring it all together
plt.plot(mtf_axis_px, mtf_optics, label='optics')
#plt.plot(mtf_axis_px[::2], mtf_optics_mo[::2], 'o', label='optics model')
plt.plot(mtf_axis_px, mtf_px, label='pixel')
#plt.plot(mtf_axis_px, mtf_px_mo,'.', label='pixel model')
plt.plot(mtf_axis_px, mtf_sys, label='system')
plt.plot(mtf_axis_px[::2], mtf_sys_mo[::2], 'o', label='pixelated model')
plt.xlabel('cycles/px')
plt.ylabel('MTF')
plt.legend(prop={'size': 12})
plt.savefig('../../_static/img/pixel_mtf.png', transparent=False, bbox_inches='tight', dpi=150)
