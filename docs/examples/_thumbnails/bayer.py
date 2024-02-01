import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

gaussian = scipy.signal.windows.gaussian(250, 25)
gaussian = gaussian[25:-25]
qe_wave = np.arange(350, 750)*1e-9

qe_red = np.zeros(qe_wave.size)
qe_green = np.zeros(qe_wave.size)
qe_blue = np.zeros(qe_wave.size)

qe_blue[10:210] = gaussian
qe_green[90:290] = gaussian
qe_red[170:370] = gaussian


plt.plot(qe_wave, qe_blue, '#5E81AC')
plt.fill_between(qe_wave, qe_blue, color='#5E81AC', alpha=0.4)

plt.plot(qe_wave, qe_green, '#A3BE8C')
plt.fill_between(qe_wave, qe_green, color='#A3BE8C', alpha=0.4)

plt.plot(qe_wave, qe_red, '#BF616A')
plt.fill_between(qe_wave, qe_red, color='#BF616A', alpha=0.4)

plt.axis('off')

plt.savefig('bayer.png', dpi=150, bbox_inches='tight')