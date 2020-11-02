from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors

POWER_STEP = 1.1 # 1 is the minimum

ALPHA_L = 1.94
L_EXCESS = 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN_RANGE=[1.0e28, 1.0e34]
L_MAX_RANGE=[1.0e34, 1.0e36]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

DRAW_EXTRA_CONTOURS = False

dimMin= int(log(L_MIN_RANGE[1]/L_MIN_RANGE[0]) / log(POWER_STEP))
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))

paperPoint = [1e35, 1e29]

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def getNumPulsars(lMin, lMax):
    lumExp = lMax**(2-ALPHA_L) * Gamma(2-ALPHA_L, lMin / lMax)
    Aexp = L_EXCESS / lumExp
    Nexp = Aexp *  lMax**(1-ALPHA_L) * Gamma(1-ALPHA_L, lMin / lMax)

    return Nexp

def getNumPulsarsAboveThreshold(lMin, lMax):
    lumExp = lMax**(2-ALPHA_L) * Gamma(2-ALPHA_L, lMin / lMax)
    Aexp = L_EXCESS / lumExp
    nAbove = Aexp *  lMax**(1-ALPHA_L) * Gamma(1-ALPHA_L, L_THRESH / lMax)
    return nAbove

def getFracLumAboveThreshold(lMin, lMax):
    # return (lum above thresh) / totalLum
    fracAbove = Gamma(2-ALPHA_L, L_THRESH / lMax) / (Gamma(2-ALPHA_L, lMin / lMax))
    return fracAbove


print("Paper values:", getNumPulsars(1e29, 1e35))
print("Delta function values:", getNumPulsars(L_THRESH, L_THRESH*1.0001))


def getGreenYIntercept():
    ERROR_WIDTH = 0.01

    lowLmin = L_MIN_RANGE[0] # Too few pulsars above threshold
    highLmin = L_MIN_RANGE[1] # Too many pulsars above threshold
    numLow = getNumPulsarsAboveThreshold(lowLmin, L_THRESH)
    numHigh = getNumPulsarsAboveThreshold(highLmin, L_THRESH)
    errorNumPulsars = 2 * ERROR_WIDTH
    while abs(numLow - numHigh) > ERROR_WIDTH:
        mid = (lowLmin + highLmin) / 2
        numMid = getNumPulsarsAboveThreshold(mid, L_THRESH)
        if numMid > NUM_PULSARS_ABOVE_THRESHOLD:
            highLmin = mid
            numHigh = numMid
        else:
            lowLmin = mid
            numLow = numMid
    
    print("Lmin for number requirement and Lmax = Lthresh: ", (lowLmin + highLmin) / 2)
    print("Number of total pulsars for number requirement and Lmax = Lthresh: ", getNumPulsars((lowLmin + highLmin) / 2, L_THRESH))

numPulsars = []
for i in range(dimMin):
    line = []
    for j in range(dimMax):
        nExp = getNumPulsars(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    numPulsars.append(line)

numAboveThreshold = []
for i in range(dimMin):
    line = []
    for j in range(dimMax):
        nExp = getNumPulsarsAboveThreshold(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    numAboveThreshold.append(line)

fracAboveThreshold = []
for i in range(dimMin):
    line = []
    for j in range(dimMax):
        nExp = getFracLumAboveThreshold(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    fracAboveThreshold.append(line)

lMinVals = [L_MIN_RANGE[0] * POWER_STEP**i for i in range(dimMin)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]


getGreenYIntercept()


fig, ax = plt.subplots(figsize=(6,4))
plt.text(0.95, 0.95, 'Greens: number limit\nReds: luminosity limit\nBold: Fermilab observations', 
    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, color='white', backgroundcolor=(0, 0, 0, 0.3))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Lmax")
plt.ylabel("Lmin")
plt.title("Exp cutoff (alpha={0})".format(ALPHA_L))

c1 = plt.pcolor(lMaxVals, lMinVals, numPulsars, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='PuBu_r')
plt.colorbar(c1, extend='max')

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, numAboveThreshold, [10*i for i in range(1, 20)], 
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)], linewidths=1)
plt.contour(lMaxVals, lMinVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[(0, 1, 0)], linewidths=2)

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, fracAboveThreshold, [0.1*i for i in range(1, 15)], 
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)], linewidths=1)
plt.contour(lMaxVals, lMinVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], colors=[(1, 0, 0)], linewidths=2)


plt.scatter(paperPoint[0], paperPoint[1], c='purple')


if(DRAW_EXTRA_CONTOURS):
    plt.savefig("contour-overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("contour-overlay.png")

plt.show()