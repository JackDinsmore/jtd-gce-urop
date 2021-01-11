''' Things to check:
    Is the NFW profile a number density?
    What is r_s from the Fermilab paper?
    I should integrate the square of the NFW profile, not square the integral, right?
    Am I using the number position distribution and luminosity position distribution in the correct places?
'''

from matplotlib import pyplot as plt
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
import numpy as np
from numpy import exp, log
from math import pi, sqrt, sin, cos
from astropy.io import fits

plt.style.use('latex')


POWER_STEP = 1.1 # 1 is the minimum

ALPHA = 1.94
L_EXCESS = 6.37e36  # All units are in ergs per second
DIST_TO_CENTER = 8.5 # kpc
CM_PER_KPC = 3.086e21
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e36]
LON_RANGE = [-20, 20]
LAT_RANGE = [-20, 20]

# Values for the computation of the integral of the NFW profile.
GAMMA = 1.2
NFW_SCALE_DIST = 1 ### FIX # kpc
INTEGRAL_HALF_WIDTH = DIST_TO_CENTER / 2 # kpc
INTEGRAL_STEPS = 10000
LOAD_DISTROS = False

SENSITIVITY_PATH = "../../sensitivity/detthresh_P8R3_source_10years_PL22.fits"

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0
DISPLAY_SIZE = 20.0

DRAW_EXTRA_CONTOURS = False

paperPoint = [1e35, 1e29]

dimMin = int(log(L_MIN_RANGE[1]/L_MIN_RANGE[0]) / log(POWER_STEP))
dimMax = int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))
plotShape = (dimMin, dimMax)

class SensitivityMap:
    def __init__(self):
        hdu_list = fits.open(SENSITIVITY_PATH)
        self.flux = hdu_list[0].data
        hdu_list.close()
        self.luminosity = self.fluxToThresholdLuminosity(self.flux)
        self.upperLeft = ( int((1 - DISPLAY_SIZE / 90) * self.flux.shape[0] / 2), int((1 - DISPLAY_SIZE / 180) * self.flux.shape[1] / 2) )
        self.lowerRight = ( int((1 + DISPLAY_SIZE / 90) * self.flux.shape[0] / 2), int((1 + DISPLAY_SIZE / 180) * self.flux.shape[1] / 2) )
        self.skyShape = (self.lowerRight[0] - self.upperLeft[0], self.lowerRight[1] - self.upperLeft[1])

    def fluxToThresholdLuminosity(self, flux):
        return flux * (4 * pi * (DIST_TO_CENTER * CM_PER_KPC) ** 2)

    def getThreshold(self, lat, lon):
        x, y = self.latLonToIndex(lat, lon)
        return self.luminosity[x + self.upperLeft[0]][y + self.upperLeft[1]]

    def indexToLatLon(self, x, y):
        deltaXFrac = 1 - x / (self.skyShape[0] / 2)
        deltaYFrac = y / (self.skyShape[1] / 2) - 1
        return DISPLAY_SIZE * deltaXFrac * pi / 180, DISPLAY_SIZE * deltaYFrac * pi / 180

    def latLonToIndex(self, lat, lon):
        deltaXFrac = lat * 180 / pi / DISPLAY_SIZE
        deltaYFrac = lon * 180 / pi / DISPLAY_SIZE
        return int((1 - deltaXFrac) * (self.skyShape[0] / 2)), int((1 + deltaYFrac) * (self.skyShape[1] / 2))

thresholds = SensitivityMap()

# ================= These functions characterize the luminosity function ============

def luminosityFunction(luminosity, lMin, lMax):
    return luminosity**(-ALPHA) * np.exp(-luminosity / lMax)
    
def integrate(start, stop, lMin, lMax): # Integrate luminosityFunction from start to stop
    return lMax ** (1 - ALPHA) * (Gamma(1 - ALPHA, start / lMax) - Gamma(1 - ALPHA, stop / lMax))

def lintegrate(start, stop, lMin, lMax): # Integrate luminosityFunction times luminosity from start to stop
    return lMax ** (2 - ALPHA) * (Gamma(2 - ALPHA, start / lMax) - Gamma(2 - ALPHA, stop / lMax))

# ================= Declare helper functions =======================

def generateNFWMaps():
    numberDistro = np.zeros(thresholds.skyShape)
    lumDistro = np.zeros(thresholds.skyShape)
    if LOAD_DISTROS:
        numberFile = open("data/nfw-number-profile.txt", 'r')
        lumFile = open("data/nfw-lum-profile.txt", 'r')
        numberText = numberFile.read()
        lumText = lumFile.read()
        numberFile.close()
        lumFile.close()
        x = 0
        for line in numberText.split(":"):
            y = 0
            for entry in line.split(", "):
                numberDistro[x][y] = float(entry)
                y += 1
            x+= 1
        x = 0
        for line in lumText.split(":"):
            y = 0
            for entry in line.split(", "):
                lumDistro[x][y] = float(entry)
                y += 1
            x+= 1

    else:
        numberText = ''
        lumText = ''
        deltaRadialDistance = INTEGRAL_HALF_WIDTH * 2 / INTEGRAL_STEPS
        for x in range(numberDistro.shape[0]):
            print(x, "/", numberDistro.shape[0])
            for y in range(numberDistro.shape[1]):
                lat, lon = thresholds.indexToLatLon(x, y)
                radialDistance = DIST_TO_CENTER - INTEGRAL_HALF_WIDTH
                numberIntegral = 0
                lumIntegral = 0
                while radialDistance < DIST_TO_CENTER + INTEGRAL_HALF_WIDTH:
                    distFromCenter = sqrt((radialDistance * cos(lon) * cos(lat) - DIST_TO_CENTER)**2 + (radialDistance * sin(lon) * cos(lat))**2 + (radialDistance * sin(lat))**2)
                    volumeElement = deltaRadialDistance * radialDistance **2
                    nfwValue = ((distFromCenter / NFW_SCALE_DIST)**-GAMMA * (1 + distFromCenter / NFW_SCALE_DIST) ** (-3 + GAMMA))**2
                    numberIntegral += nfwValue * deltaRadialDistance * radialDistance**2 # * deltaLon * deltaLat #Ignore the angular part since it introduces a constant cofactor.
                    lumIntegral += nfwValue * deltaRadialDistance * DIST_TO_CENTER**2 # * deltaLon * deltaLat #Ignore the angular part since it introduces a constant cofactor.
                    radialDistance += deltaRadialDistance
                numberDistro[x][y] = numberIntegral
                lumDistro[x][y] = lumIntegral
                numberText += str(numberIntegral)
                lumText += str(lumIntegral)
                if y != numberDistro.shape[1] - 1:
                    numberText += ", "
                    lumText += ", "
            if x != numberDistro.shape[0] - 1:
                numberText += ":"
                lumText += ":"
        
        numberFile = open("data/nfw-number-profile.txt", 'w')
        lumFile = open("data/nfw-lum-profile.txt", 'w')
        numberFile.write(numberText)
        lumFile.write(lumText)
        numberFile.close()
        lumFile.close()

    return numberDistro, lumDistro


def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def totalNumFunc(threshold, lMin, lMax): # Returns (unscaled) number of pulsars
    return integrate(lMin, lMax, lMin, lMax)

def numSeenFunc(threshold, lMin, lMax): # Returns (unscaled) number of visible pulsars
    assert(threshold > lMin)
    return integrate(threshold, lMax, lMin, lMax)

def lumSeenFunc(threshold, lMin, lMax): # Returns (unscaled) fraction of luminosity visible
    assert(threshold > lMin)
    return lintegrate(threshold, lMax, lMin, lMax)

def totalLumFunc(threshold, lMin, lMax): # Returns (unscaled) total luminosity
    return lintegrate(lMin, lMax, lMin, lMax)

def getValueAtPos(lat, lon, valueFunc, lMin, lMax): # returns the (unscaled) value for a specific lon and lat 
    threshold = thresholds.getThreshold(lat, lon)
    return valueFunc(threshold, lMin, lMax)

def generateValueSkyMap(valueFunc, positionDistro, lMin, lMax): # returns an array of numTotal, numSeen, fracSeen across the entire sky
    # The values are scaled relative to each other in the image, but do not produce the correct total luminosity across the entire GCE
    skyMap = np.zeros(thresholds.skyShape)

    for x in range(skyMap.shape[0]):
        for y in range(skyMap.shape[1]):
            lat, lon = thresholds.indexToLatLon(x, y)
            val = positionDistro[x][y] * getValueAtPos(lat, lon, valueFunc, lMin, lMax)
            skyMap[x][y] = val

    return skyMap

def getValueAtConfig(valueFunc, positionDistro, lMin, lMax): # Return the (unscaled) value at a certain lMin or lMax
    skyMap = generateValueSkyMap(valueFunc, positionDistro, lMin, lMax)
    return np.sum(skyMap)

def generatePlotMap(valueFunc, positionDistro): # Return the (unscaled) plot map
    dataMap = np.zeros(plotShape)
    for i in range(dimMin):
        lMin = L_MIN_RANGE[0] * POWER_STEP**i
        print(i, '/', dimMin)
        for j in range(dimMax):
            lMax = L_MAX_RANGE[0] * POWER_STEP**i
            dataMap[i][j] = getValueAtConfig(valueFunc, positionDistro, lMin, lMax)
    return dataMap

# ========================= Generate data =========================
numberDistro, lumDistro = generateNFWMaps()
print("Distros generated")
paperScale = L_EXCESS / getValueAtConfig(totalLumFunc, lumDistro, paperPoint[1], paperPoint[0])
print("""Paper values: (lMin = {0} ergs/s, lMax = {1} ergs/s)
    Total number of pulsars: {2}
    Number of visible pulsars: {3}
    Fraction of seen luminosity: {4}
    Total luminosity: {5} ergs/s""".format(paperPoint[1], paperPoint[0],
        getValueAtConfig(totalNumFunc, numberDistro, paperPoint[1], paperPoint[0]) * paperScale,
        getValueAtConfig(numSeenFunc, numberDistro, paperPoint[1], paperPoint[0]) * paperScale,
        getValueAtConfig(lumSeenFunc, lumDistro, paperPoint[1], paperPoint[0]) * paperScale / L_EXCESS,
        getValueAtConfig(totalLumFunc, lumDistro, paperPoint[1], paperPoint[0]) * paperScale,
    ))

# Generate unscaled data
totalNum = generatePlotMap(totalNumFunc, numberDistro)
numSeen = generatePlotMap(numSeenFunc, numberDistro)
lumSeen = generatePlotMap(lumSeenFunc, lumDistro)
totalLum = generatePlotMap(totalLumFunc, lumDistro)

# Scale the plots
scale = L_EXCESS / totalLum
totalNum *= scale
numSeen *= scale
lumSeen *= scale

 # ========================= Display data =========================

lMinVals = [L_MIN_RANGE[0] * POWER_STEP**i for i in range(dimMin)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]

fig, ax = plt.subplots(figsize=(6,4))
plt.text(0.95, 0.95, 'Green: number limit\nRed: luminosity limit', 
    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, color='white', backgroundcolor=(0, 0, 0, 0.3))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Lmax")
plt.ylabel("Lmin")
plt.title("Exp cutoff (alpha={0})".format(ALPHA))

c1 = plt.pcolor(lMaxVals, lMinVals, totalNum, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='PuBu_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("Number of MSPs")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, numSeen, [10*i for i in range(1, 20)], 
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)], linewidths=1)
plt.contour(lMaxVals, lMinVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[(0, 1, 0)], linewidths=2)

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, lumSeen, [0.1*i for i in range(1, 15)], 
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)], linewidths=1)
plt.contour(lMaxVals, lMinVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[(1, 0, 0)], linewidths=2)


plt.scatter(paperPoint[0], paperPoint[1], c='purple')


if(DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay.png")

plt.show()