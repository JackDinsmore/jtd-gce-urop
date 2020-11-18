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
hdu_list.info()

image_data = hdu_list[0].data
print(hdu_list[0].header)

print(type(image_data))
print(image_data.shape)

plt.figure(figsize=(10,5))
c = plt.imshow(image_data, 
                   norm=colors.LogNorm(vmin=np.percentile(image_data, 7),
                   vmax=np.percentile(image_data, 99)))
plt.colorbar(c)
plt.title("Fermi LAT sensitivity")
plt.savefig("sensitivity.png")
plt.show()


hdu_list.close()