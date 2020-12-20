from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from math import log10

plt.style.use('latex')

hdul = fits.open('gll_psc_v27.fit')
#hdul.info()
#print('\n'.join(str(hdul[1].header).split("/")))
data = hdul[1].data
hdul.close()

DIST_FROM_GALACTIC_CENTER_KILOPARSEC = 8.5
CM_PER_KILOPARSEC = 3.086e21
PI = 3.1415926535



excludeMask = np.ones(data.shape[0], dtype=bool)


# Functions to filter point sources I don't want to include in my analysis
def excludeLatLong(excludeMask):
    thisMask = np.ones(data.shape[0], dtype=bool)
    thisMask = thisMask & ((data["GLON"] < 10) | (data["GLON"] > 350))
    thisMask = thisMask & (-10 < data["GLAT"]) & (data["GLAT"] < 10)
    return excludeMask & thisMask

def excludeDistance(excludeMask):
    thisMask = np.ones(data.shape[0], dtype=bool)
    return excludeMask & thisMask

def convertFluxToLuminosity(flux):
    return flux * (4 * PI * (DIST_FROM_GALACTIC_CENTER_KILOPARSEC * CM_PER_KILOPARSEC)**2)


# Apply filters
excludeMask = excludeLatLong(excludeMask)
excludeMask = excludeDistance(excludeMask)

print(data.shape)
data = data[excludeMask]
print(data.shape)


# Generate flux histogram
fluxData = data["Energy_Flux100"]
bins = 10**np.linspace(start=log10(np.min(fluxData)), stop=log10(np.max(fluxData)), num=50)

plt.figure()
plt.hist(fluxData, bins=bins)
plt.xlabel("Flux received (ergs / cm^2 / s)")
plt.ylabel("Point source count")
plt.xscale('log')
plt.title("Histogram of point source fluxes (0.1-100 GeV range)")


# Generate luminosity histogram
lumData = convertFluxToLuminosity(fluxData)
bins = 10**np.linspace(start=log10(np.min(lumData)), stop=log10(np.max(lumData)), num=50)

plt.figure()
fluxData = data["Energy_Flux100"]
plt.hist(lumData, bins=bins)
plt.xlabel("Luminosity (ergs / s)")
plt.ylabel("Point source count")
plt.xscale('log')
plt.axvline(x=1e34, ymin=0, ymax=data.shape[0], color='k')
plt.title("Histogram of point source luminosities (0.1-100 GeV range)")
plt.show()