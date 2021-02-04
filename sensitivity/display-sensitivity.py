# Code extracted from this tutorial:
# https://learn.astropy.org/FITS-images.html

# Data obtained from here:
# https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.utils.data import download_file
from astropy.io import fits

image_file = "detthresh_P8R3_source_10years_PL22.fits"

hdu_list = fits.open(image_file)
#hdu_list.info()

image_data = hdu_list[0].data
print('\n'.join(str(hdu_list[0].header).split('/ ')))

print(type(image_data))
print(image_data.shape)

plt.figure(figsize=(10,5))
c = plt.imshow(image_data, 
                   norm=colors.LogNorm(vmin=np.percentile(image_data, 7),
                   vmax=np.percentile(image_data, 99)))
cbar = plt.colorbar(c)
cbar.set_label("erg/cm^2/s")

plt.title("Fermi LAT sensitivity (0.1-100 GeV)")
plt.savefig("sensitivity.png")
NUM_X_LABELS = 8
xPositions = np.arange(0,image_data.shape[1],image_data.shape[1]//NUM_X_LABELS) # pixel count at label position
xLabels = np.linspace(start=-180, stop=180, num=NUM_X_LABELS+1)[:-1] # labels you want to see
xPositions = xPositions[:NUM_X_LABELS]
plt.xticks(xPositions, xLabels)

NUM_Y_LABELS = 5
yPositions = np.arange(0,image_data.shape[0],image_data.shape[0]//NUM_Y_LABELS) # pixel count at label position
yLabels = np.linspace(start=-90, stop=90, num=NUM_Y_LABELS+1)[:-1] # labels you want to see
yPositions = yPositions[:NUM_Y_LABELS]
plt.yticks(yPositions, yLabels)

plt.xlabel("latitude (deg)")
plt.ylabel("longitude (deg)")
plt.show()


hdu_list.close()