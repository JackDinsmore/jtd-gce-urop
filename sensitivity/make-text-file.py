# FITS file code extracted from this tutorial:
# https://learn.astropy.org/FITS-images.html

# Flux data obtained from here:
# https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/

# NFW profile formula obtained from here:
# https://arxiv.org/pdf/1911.12369.pdf

import numpy as np
from astropy.utils.data import download_file
from astropy.io import fits

image_file = "detthresh_P8R3_source_10years_PL22.fits"

hdu_list = fits.open(image_file)
#hdu_list.info()
image_data = hdu_list[0].data
hdu_list.close()

f = open("sensitivity.txt", 'w')
x = 0
for line in image_data:
    f.write(','.join([str(item) for item in line]))
    if x != image_data.shape[0] - 1:
        f.write("\n")
    x+= 1
f.close()