import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits

NIRCAM_FITS = '/Users/akee/Downloads/webbpsf-data/NIRCam/OPD/OPD_RevW_ote_for_NIRCam_predicted.fits'

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

if NIRCAM_FITS:
    hdul = fits.open(NIRCAM_FITS)
    fig, ax = plt.subplots()
    ax.imshow(hdul[0].data[0], origin='lower')

    fig.savefig('../../../_static/img/nircam.png', dpi=dpi*2)
