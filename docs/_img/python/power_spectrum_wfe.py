import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.circle((256, 256), 128)
w = lentil.power_spectrum(amp, 1/100, 5e-9, 8, 3)
plt.imshow(w, origin='lower')
plt.colorbar()
plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.savefig('../../_static/img/power_spectrum_wfe.png', transparent=True, bbox_inches='tight', dpi=150)
