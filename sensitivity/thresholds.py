# FITS file code extracted from this tutorial:
# https://learn.astropy.org/FITS-images.html

# Flux data obtained from here:
# https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

# NFW profile formula obtained from here:
# https://arxiv.org/pdf/1911.12369.pdf

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.utils.data import download_file
from astropy.io import fits
from math import pi, cos, sqrt

matplotlib.font_manager.FontManager()
plt.style.use('latex')

DIST_TO_CENTER = 8.5
KAPPA = 0.5 # dist to center / r_s
GAMMA = 1.2 # Shape facter of the NFW profile
A = 1 # Something like the number of objects per cubic kiloparsec; in fact, it's the coefficient of the NFW number density
CM_PER_KPC = 3.086e21

DISPLAY_SIZE = 20.0
THRESHOLD = 1e34

def fluxToThresholdLuminosity(flux):
    return flux * (4 * pi * (DIST_TO_CENTER * CM_PER_KPC) ** 2)

image_file = "detthresh_P8R3_source_10years_PL22.fits"

hdu_list = fits.open(image_file)
#hdu_list.info()
image_data = hdu_list[0].data
hdu_list.close()

lum_data = np.zeros_like(image_data)
print(image_data.shape)
deltaLat = pi / image_data.shape[0]
deltaLon = 2 * pi / image_data.shape[1]
print(deltaLat, deltaLon)

trimmed_flux_data = []
trimmed_lum_data = []
for x in range(len(image_data)):
    #print("Latitude:", x * deltaLat * 180 / pi - 90)
    fluxLine = []
    lumLine = []
    for y in range(len(image_data[x])):
        flux = image_data[x][y]
        lat = x * deltaLat - pi / 2
        lon = y * deltaLon - pi
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
print(trimmed_flux_data)
print(trimmed_lum_data)


# Display flux data:
plt.figure()
c = plt.imshow(trimmed_flux_data, 
                   norm=colors.LogNorm(vmin=np.percentile(trimmed_flux_data, 7),
                   vmax=np.percentile(trimmed_flux_data, 99)))
cbar = plt.colorbar(c)
cbar.set_label("erg/cm^2/s")

plt.title("Fermi LAT flux sensitivity (0.1-100 GeV)")
NUM_X_LABELS = 8
xPositions = np.arange(0, trimmed_flux_data.shape[1], trimmed_flux_data.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_X_LABELS+1)[:-1] # labels you want to see
xPositions = xPositions[:NUM_X_LABELS]
plt.xticks(xPositions, xLabels)

NUM_Y_LABELS = 5
yPositions = np.arange(0, trimmed_flux_data.shape[0], trimmed_flux_data.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_Y_LABELS+1)[:-1] # labels you want to see
yPositions = yPositions[:NUM_Y_LABELS]
plt.yticks(yPositions, yLabels)
plt.xlabel("latitude (deg)")
plt.ylabel("longitude (deg)")
plt.tight_layout()
plt.savefig("sensitivity.png")

# Cut luminosity data
'''for x in range(len(trimmed_lum_data)):
    for y in range(len(trimmed_lum_data[x])):
        if trimmed_lum_data[x][y] > THRESHOLD:
            trimmed_lum_data[x][y] = THRESHOLD'''


# Display luminosity data

plt.figure()
c = plt.imshow(trimmed_lum_data, 
                   norm=colors.LogNorm(vmin=np.percentile(trimmed_lum_data, 7),
                   vmax=np.percentile(trimmed_lum_data, 99)))
cbar = plt.colorbar(c)
cbar.set_label("erg/s")

plt.title("Fermi LAT luminosity threshold (0.1-100 GeV)")
xPositions = np.arange(0, trimmed_lum_data.shape[1], trimmed_lum_data.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_X_LABELS+1)[:-1] # labels you want to see
xPositions = xPositions[:NUM_X_LABELS]
plt.xticks(xPositions, xLabels)
yPositions = np.arange(0, trimmed_lum_data.shape[0], trimmed_lum_data.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=-DISPLAY_SIZE, stop=DISPLAY_SIZE, num=NUM_Y_LABELS+1)[:-1] # labels you want to see
yPositions = yPositions[:NUM_Y_LABELS]
plt.yticks(yPositions, yLabels)

plt.xlabel("latitude (deg)")
plt.ylabel("longitude (deg)")
plt.tight_layout()
plt.savefig("luminosity-thresholds.png")
plt.show()

