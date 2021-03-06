from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('latex')

POWER_STEP = 1.1 # 1 is the minimum

FERMILAB_GNFW_FLUX = 5.716826056794996e-10
CALORE_FLUX = 1.3592975571782718e-09
DI_MAURO_FLUX = 1.2953417255755896e-09
AJELLO_FLUX = 1.8268259089431682e-09
FERMILAB_NFW_FLUX = 7.770902203635269e-10
ABAZAJIAN_FLUX = 4.887121798902831e-10

ALPHA_L = 1.94
L_EXCESS = 1.2953417255755896e-09 / 8.331593765023139e-47#6.756e36# 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

DRAW_EXTRA_CONTOURS = True
SHOW_NUMBERS = False
LINE_COLOR = (0.8, 0.3, 0.1)
SHADE_SCALE=25

dimMin= int(log(L_MIN_RANGE[1]/L_MIN_RANGE[0]) / log(POWER_STEP))
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))

paperPoint = [1e35, 1e29]

def shade(field, threshold, xs, ys, off=False):
    px = []
    py = []
    for x in range(0 if off else 1, SHADE_SCALE, 1):
        inx = int(float(x) / SHADE_SCALE * field.shape[1])
        for y in range(0 if off else 1, SHADE_SCALE, 1):
            iny = int(float(y) / SHADE_SCALE * field.shape[0])
            if field[iny][inx] < threshold:
                fracx = float(x) / SHADE_SCALE * field.shape[1] - inx
                fracy = float(y) / SHADE_SCALE * field.shape[0] - iny
                px.append(xs[inx] + fracx * (xs[inx+1] - xs[inx]))
                py.append(ys[iny] + fracy * (ys[iny+1] - ys[iny]))
    plt.scatter(px, py, marker=('|' if off else '_'), c='black', sizes = (20,), alpha=0.5)
        

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


print("""Paper point info:
    Coordinates: L0 = {0}, sigma = {1}
    Total num pulsars: {2}
    Num pulsars above threshold: {3}
    Fraction luminosity above threshold: {4}""".format(paperPoint[1], paperPoint[0], getNumPulsars(paperPoint[1], paperPoint[0]), 
    getNumPulsarsAboveThreshold(paperPoint[1], paperPoint[0]), getFracLumAboveThreshold(paperPoint[1], paperPoint[0])))

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
    numAboveThreshold.append(np.asarray(line))

fracAboveThreshold = []
for i in range(dimMin):
    line = []
    for j in range(dimMax):
        nExp = getFracLumAboveThreshold(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    fracAboveThreshold.append(np.asarray(line))

numAboveThreshold = np.stack(numAboveThreshold)
fracAboveThreshold = np.stack(fracAboveThreshold)



lMinVals = [L_MIN_RANGE[0] * POWER_STEP**i for i in range(dimMin)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]


getGreenYIntercept()


fig, ax = plt.subplots(figsize=(6, 4))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("$L_{max}$")
plt.ylabel("$L_{min}$")
plt.title("Power law, step function")


if SHOW_NUMBERS:
    c1 = plt.pcolor(lMaxVals, lMinVals, numPulsars, 
                    norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                    vmax=max([max(v) for v in numPulsars])), cmap='Greys_r')
    cbar = plt.colorbar(c1, extend='max')
    cbar.set_label("$N$")


# Greens
if DRAW_EXTRA_CONTOURS:
    plt.contour(lMaxVals, lMinVals, numAboveThreshold * FERMILAB_GNFW_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["blue"], linewidths=1, label="$N_r$ Fermilab gNFW")
    plt.contour(lMaxVals, lMinVals, numAboveThreshold * ABAZAJIAN_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["green"], linewidths=1, label="$N_r$ Abazajian")
    plt.contour(lMaxVals, lMinVals, numAboveThreshold * AJELLO_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["teal"], linewidths=1, label="$N_r$ Ajello")
plt.contour(lMaxVals, lMinVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR], linewidths=2, label="$N_r=47$")

# Reds
plt.contour(lMaxVals, lMinVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2, label="$R_r=0.2$")

# Observation
if DRAW_EXTRA_CONTOURS:
    shade(numAboveThreshold * AJELLO_FLUX / DI_MAURO_FLUX, NUM_PULSARS_ABOVE_THRESHOLD, lMaxVals, lMinVals)
else:
    shade(numAboveThreshold, NUM_PULSARS_ABOVE_THRESHOLD, lMaxVals, lMinVals)
shade(fracAboveThreshold, FRAC_ABOVE_THRESHOLD, lMaxVals, lMinVals, True)


# Final points 
plt.scatter(paperPoint[0], paperPoint[1], c='purple')

custom_lines = [Line2D([0], [0], color=LINE_COLOR, lw=2),
                Line2D([0], [0], color="blue", lw=1),
                Line2D([0], [0], color="green", lw=1),
                Line2D([0], [0], color="teal", lw=1),
                Line2D([0], [0], color=LINE_COLOR, lw=2, dashes=(4, 2)),]
plt.legend(custom_lines, ['$N_r$ Di Mauro', "$N_r$ Fermilab gNFW", "$N_r$ Abazajian", "$N_r$ Ajello", '$R_r$'], loc="upper right")
plt.xlim(lMaxVals[0], lMaxVals[-1])
plt.ylim(lMinVals[0], lMinVals[-1])
plt.tight_layout()

print([np.percentile(numPulsars, i) for i in range(0, 100, 10)])


if(DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay.png")

plt.show()