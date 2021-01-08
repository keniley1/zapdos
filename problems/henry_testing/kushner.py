import numpy as np
import matplotlib.pyplot as plt


D_g = 1.0
D_l = 1e-3
dx = 1.0


# Henry coefficient (for OH in this case)
h = 6.2e2

n_g = 1
#n_l = np.linspace(1e-4, 1e4, 1000)
#n_l = np.logspace(-4, 2.8, 1000)
n_l = np.logspace(0, 2.8, 1000)

diff_gl = n_g - n_l
henry_gl = (1.0 - n_l/(n_g * h))

diff_lg = n_l - n_g
henry_lg = (1.0 - n_g/(n_l / h))

#flux = D_g / dx * (1.0 - n_l/(n_g * h)) * (n_g - n_l)
flux_gl = D_g / dx * henry_gl * diff_gl
flux_lg = D_l / dx * henry_lg * diff_lg 

for i in range(len(n_l)):
    if (n_l[i] <= (n_g * h)):
        flux_lg[i] = 0
    else:
        flux_gl[i] = 0

#plt.semilogx(n_l, flux)
fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(8,12))
ax1.semilogx(n_l, flux_gl)
ax1.semilogx(n_l, flux_lg)

ax2.semilogx(n_l, henry_gl)

ax3.semilogx(n_l, diff_gl)

plt.show()

