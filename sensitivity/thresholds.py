# FITS file code extracted from this tutorial:
# https://learn.astropy.org/FITS-images.html

# Flux data obtained from here:
# https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

# NFW profile formula obtained from here:
# https://arxiv.org/pdf/1911.12369.pdf

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.utils.data import download_file
from astropy.io import fits
from math import pi, cos, sqrt

plt.style.use('jcap')

DIST_TO_CENTER = 8.5
KAPPA = 0.5 # dist to center / r_s
GAMMA = 1.2 # Shape facter of the NFW profile
A = 1 # Something like the number of objects per cubic kiloparsec; in fact, it's the coefficient of the NFW number density
CM_PER_KPC = 3.086e21
FLUX_TO_LUM = 1.1093417307914119e-46

DISPLAY_SIZE = 20.0
THRESHOLD = 1e34
NUM_X_LABELS = 8
NUM_Y_LABELS = 8

def fluxToThresholdLuminosity(flux):
    return flux / FLUX_TO_LUM #flux * (4 * pi * (DIST_TO_CENTER * CM_PER_KPC) ** 2)

image_file = "detthresh_P8R3_source_10years_PL22.fits"

hdu_list = fits.open(image_file)
#hdu_list.info()
image_data = hdu_list[0].data
hdu_list.close()

lum_data = np.zeros_like(image_data)
deltaLat = pi / image_data.shape[0]
deltaLon = 2 * pi / image_data.shape[1]

trimmed_flux_data = []
trimmed_lum_data = []
for x in range(len(image_data)):
    fluxLine = []
    lumLine = []
    for y in range(len(image_data[x])):
        flux = image_data[x][y]
        lat = x * deltaLat - pi / 2
        lon = y * deltaLon - pi
        if abs(lat) < 2 * pi / 180:
            flux = np.nan#1e-12
        if abs(lat) > DISPLAY_SIZE * pi / 180 or abs(lon) > DISPLAY_SIZE * pi / 180:
            lum_data[x][y] = 0.1
        else:
            lum = fluxToThresholdLuminosity(flux)
            lum_data[x][y] = lum
            fluxLine.append(flux)
            lumLine.append(lum)
    if len(fluxLine) > 0:
        trimmed_flux_data.append(np.asarray(fluxLine))
        trimmed_lum_data.append(np.asarray(lumLine))

trimmed_flux_data = np.asarray(trimmed_flux_data)
trimmed_lum_data = np.asarray(trimmed_lum_data)

#plt.figure(figsize=(8, 6))
#mpl.rcParams["font.size"]=12

# Display flux data:
c = plt.imshow(trimmed_flux_data,
                   norm=colors.LogNorm(vmin=np.nanmin(trimmed_flux_data),
                   vmax=np.nanpercentile(trimmed_flux_data, 99)))
cbar = plt.colorbar(c)
cbar.set_label("$F_\\textrm{th}(b, \\ell)$ (erg/cm$^2$/s)")

xPositions = np.arange(0, trimmed_flux_data.shape[1], trimmed_flux_data.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_X_LABELS+1) # labels you want to see
np.append(xPositions, trimmed_flux_data.shape[1])
plt.xticks(xPositions, [int(i) for i in xLabels])

yPositions = np.arange(0, trimmed_flux_data.shape[0], trimmed_flux_data.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_Y_LABELS+1) # labels you want to see
np.append(yPositions, trimmed_flux_data.shape[0])
plt.yticks(yPositions, [int(i) for i in yLabels])
plt.xlabel("$\\ell$ (deg)")
plt.ylabel("$b$ (deg)")
plt.tight_layout()
plt.savefig("flux-thresholds.pdf")

# Cut luminosity data
'''for x in range(len(trimmed_lum_data)):
    for y in range(len(trimmed_lum_data[x])):
        if trimmed_lum_data[x][y] > THRESHOLD:
            trimmed_lum_data[x][y] = THRESHOLD'''


# Display luminosity data

plt.figure()
c = plt.imshow(trimmed_lum_data,
                   norm=colors.LogNorm(vmin=np.nanmin(trimmed_lum_data),
                   vmax=np.nanpercentile(trimmed_lum_data, 99)))
cbar = plt.colorbar(c)
cbar.set_label("$L_{th}(b, \\ell)$ (erg/s)")

plt.title("Fermi LAT luminosity threshold (0.1-100 GeV)")
xPositions = np.arange(0, trimmed_lum_data.shape[1], trimmed_lum_data.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_X_LABELS+1) # labels you want to see
np.append(xPositions, trimmed_flux_data.shape[1])
plt.xticks(xPositions, [int(i) for i in xLabels])
yPositions = np.arange(0, trimmed_lum_data.shape[0], trimmed_lum_data.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_Y_LABELS+1) # labels you want to see
np.append(yPositions, trimmed_flux_data.shape[0])
plt.yticks(yPositions, [int(i) for i in yLabels])

print(np.nanmin(trimmed_lum_data))

plt.xlabel("$\\ell$ (deg)")
plt.ylabel("$b$ (deg)")
plt.tight_layout()
plt.savefig("luminosity-thresholds.png")
plt.show()
