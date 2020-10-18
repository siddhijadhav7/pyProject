import numpy as np
import healpy as hp
from astropy import constants as const
import matplotlib.pyplot as plt

def T_o(n_o):
  T_o = T_e/(gamma*(1-np.dot(n_o,beta)))
  return T_o

v = 369.82 #km/s
l = np.radians(263.99) # Longitude
b = np.radians(48.26)  # Latitude
c = const.c.value

T_e = 2.7255 # K
beta = [v*np.sin(b)*np.cos(l),v*np.sin(b)*np.sin(l),v*np.cos(b)]
gamma = np.sqrt(1 - pow(v/c,2))


NSIDE = 32
npix = hp.nside2npix(NSIDE)
pix_centers = hp.pix2vec(NSIDE, np.arange(npix))
T = np.zeros(npix)


for i in range(npix):
  n_o = [pix_centers[0][i],pix_centers[1][i],pix_centers[2][i]]
  T[i] = T_o(n_o)

hp.mollview(T,norm="hist", title="CMB Dipole Anisotropy")
plt.show()
