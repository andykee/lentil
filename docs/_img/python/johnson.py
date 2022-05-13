import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits

# https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/johnson_*.fits

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

fig, ax = plt.subplots()

for f in ('U','B','V','R','I'):
    hdul = fits.open(f'johnson_{f.lower()}.fits')
    ax.plot(hdul[1].data['WAVELENGTH']/10, hdul[1].data['THROUGHPUT'], label=f)
    ax.grid()
    ax.legend()
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Transmission [A.U.]')

    fig.savefig('../../_static/img/johnson.png', dpi=dpi*2)
