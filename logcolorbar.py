import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colorbar import make_axes
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
N = 100
X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

fig, ax = plt.subplots(2, 1)
pcm = ax[0].pcolormesh(X, Y, Z,
                       norm=colors.SymLogNorm(linthresh=0.1, linscale=1),
                       cmap='summer')
# fig.colorbar(pcm, ax=ax[0], extend='both')
fig.colorbar(pcm, ax=ax[0], extend='both')
# pcm = ax[1].pcolormesh(X, Y, Z, cmap='RdBu_r', vmin=-np.max(Z))
# fig.colorbar(pcm, ax=ax[1], extend='both')
plt.show()