import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits

# https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/alpha_lyr_stis_010.fits
VEGA = 'alpha_lyr_stis_010.fits'

dpi = 144
mpl.rcParams['figure.figsize'] = (2.5, 2.5)
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['font.size'] = 12*72/dpi
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.transparent'] = True
mpl.rcParams['axes.titlesize'] = 14*72/dpi  # 14 pt
mpl.rcParams['axes.labelsize'] = 12*72/dpi  # 12 pt
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.major.size'] = 2
mpl.rcParams['grid.linewidth'] = 0.5
mpl.rcParams['xtick.labelsize'] = 12*72/dpi  # 12 pt
mpl.rcParams['ytick.labelsize'] = 12*72/dpi  # 12 pt
mpl.rcParams['legend.fontsize'] = 12*72/dpi  # 12 pt
mpl.rcParams['lines.linewidth'] = 1

if VEGA:
    hdul = fits.open(VEGA)
    fig, ax = plt.subplots()
    ax.plot(hdul[1].data['WAVELENGTH'][0:3680]/10, hdul[1].data['FLUX'][0:3680])
    ax.grid()
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Flux [erg s^-1 cm^-2]')

    fig.savefig('../../_static/img/vega.png', dpi=dpi*2)
