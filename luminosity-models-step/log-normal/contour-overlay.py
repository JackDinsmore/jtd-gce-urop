from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')

DIM_TRIALS=100

L_EXCESS = 1.2953417255755896e-09 / 8.331593765023139e-47# 6.756e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_0_RANGE=[1.0e30, 2.0e36]#[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]
powerStep =(L_0_RANGE[1] / L_0_RANGE[0])**(1/DIM_TRIALS)

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0
LINE_COLOR = "C1"
DRAW_EXTRA_CONTOURS =False
DRAW_PLOEG_POINT = True

paperPoint = [0.88e34, 0.62]
ploegPoint = [10**32.206, 0.70585]
SHADE_SCALE=25

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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.7)

def getTotalLum(L0, sigma):
    return exp(0.5 * sigma**2 * log(10)**2) * L0

def getNumPulsars(L0, sigma):
    scale = L_EXCESS / getTotalLum(L0, sigma)
    return scale * 1 # scale times number of pulsars per unit

def getNumPulsarsAboveThreshold(L0, sigma):
    scale = L_EXCESS / getTotalLum(L0, sigma)
    return scale * 0.5 * erfc((log(L_THRESH) - log(L0)) / (sqrt(2) * sigma * log(10)))

def getFracLumAboveThreshold(L0, sigma):
    # return (lum above thresh) / (total lum)
    erfArgument = sqrt(2) * (sigma**2 * log(10)**2 - log(L_THRESH) + log(L0)) / (sigma * log(100))
    amountAbove = 0.5 * (1 + erf(erfArgument))
    return amountAbove # See the mathematica notebook for a derivation


def getMinNumPulsarsInTriangle():
    EXTRA_DENSITY = 2
    minSoFar = -1
    powerStepNew =(L_0_RANGE[1] / L_0_RANGE[0])**(1/(DIM_TRIALS * EXTRA_DENSITY))
    for j in range(DIM_TRIALS * EXTRA_DENSITY):
        sigma= SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / (DIM_TRIALS * EXTRA_DENSITY) * j
        for i in range(DIM_TRIALS * EXTRA_DENSITY):
            L0=L_0_RANGE[0] * powerStepNew**i
            nAbove = getNumPulsarsAboveThreshold(L0, sigma)
            if nAbove > NUM_PULSARS_ABOVE_THRESHOLD:
                continue
            fracAbove = getFracLumAboveThreshold(L0, sigma)
            if fracAbove > FRAC_ABOVE_THRESHOLD:
                continue
            n = getNumPulsars(L0, sigma)
            if(n < minSoFar or minSoFar < 0):
                minSoFar = n

    print("Minimum number of total pulsars in allowed triangle: ", minSoFar)
    return minSoFar

def getMinPulsarsWithinOneStdevOfSigma():
    SHORT_SIGMA_RANGE = [0.62 - 0.16, 0.62 + 0.15]
    EXTRA_DENSITY = 2
    minPulsars = -1
    minL0 = 0
    minSigma = 0
    powerStepNew =(L_0_RANGE[1] / L_0_RANGE[0])**(1/(DIM_TRIALS * EXTRA_DENSITY))
    for j in range(DIM_TRIALS * EXTRA_DENSITY):
        sigma= SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / (DIM_TRIALS * EXTRA_DENSITY) * j
        if sigma < SHORT_SIGMA_RANGE[0] or SHORT_SIGMA_RANGE[1] < sigma:
            continue
        tooFewPulsarsAboveThreshold = False
        for i in range(DIM_TRIALS * EXTRA_DENSITY):
            L0=L_0_RANGE[0] * powerStepNew**i

            if sigma > 0.7 and L0 > 2e33:
                continue # Cut out top green line

            tooFewPulsarsAboveThresholdNow = getNumPulsarsAboveThreshold(L0, sigma) < NUM_PULSARS_ABOVE_THRESHOLD
            if tooFewPulsarsAboveThresholdNow ==(not tooFewPulsarsAboveThreshold) and i > 0:
                # We have just crossed over the green line
                numPulsarsNow = getNumPulsars(L0, sigma)
                if numPulsarsNow < minPulsars or minPulsars < 0:
                    minPulsars = numPulsarsNow
                    minL0 = L0
                    minSigma = sigma
            tooFewPulsarsAboveThreshold = tooFewPulsarsAboveThresholdNow

    print("Fewest possible pulsars required to hit the green line with sigma in 1 stdev of paper values: {0} at coordinates L_0={1}, sigma={2}".format(minPulsars, minL0, minSigma))
    return (minL0, minSigma)

def getPaperPointInfo(L0 = paperPoint[0], sigma = paperPoint[1]):
    print("""Paper point info:
    Coordinates: L0 = {0}, sigma = {1}
    Total num pulsars: {2}
    Num pulsars above threshold: {3}
    Fraction luminosity above threshold: {4}""".format(L0, sigma, getNumPulsars(L0, sigma),
    getNumPulsarsAboveThreshold(L0, sigma), getFracLumAboveThreshold(L0, sigma)))

def getPloegPointInfo(L0, sigma):
    print("""Ploeg point info:
    Coordinates: L0 = {0}, sigma = {1}
    Total num pulsars: {2}
    Num pulsars above threshold: {3}
    Fraction luminosity above threshold: {4}""".format(L0, sigma, getNumPulsars(L0, sigma),
    getNumPulsarsAboveThreshold(L0, sigma), getFracLumAboveThreshold(L0, sigma)))


getMinNumPulsarsInTriangle()
minPoint = getMinPulsarsWithinOneStdevOfSigma()
getPaperPointInfo()
getPloegPointInfo(ploegPoint[0], ploegPoint[1])


numPulsars = []
numAboveThreshold = []
fracAboveThreshold = []
for j in range(DIM_TRIALS):
    lineNumPulsars = []
    lineNumAboveThreshold = []
    lineFracAboveThreshold = []
    for i in range(DIM_TRIALS):
        L0 = L_0_RANGE[0] * powerStep**i
        sigma =  SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j

        numNow = getNumPulsars(L0, sigma)
        numAbove = getNumPulsarsAboveThreshold(L0, sigma)
        fracAbove = getFracLumAboveThreshold(L0, sigma)

        lineNumPulsars.append(numNow)
        lineNumAboveThreshold.append(numAbove)
        lineFracAboveThreshold.append(fracAbove)

    numPulsars.append(lineNumPulsars)
    numAboveThreshold.append(np.asarray(lineNumAboveThreshold))
    fracAboveThreshold.append(np.asarray(lineFracAboveThreshold))

numAboveThreshold = np.stack(numAboveThreshold)
fracAboveThreshold = np.stack(fracAboveThreshold)


xVals = [L_0_RANGE[0] * powerStep**i for i in range(DIM_TRIALS)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j for j in range(DIM_TRIALS)]


fig, ax = plt.subplots()
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])

plt.xscale("log")
plt.ylabel("$\sigma$")
plt.xlabel("$L_0$ [ergs/s]")

cols = colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars]))
c1 = plt.contourf(xVals, yVals, numPulsars,
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N_\\textrm{GCE}$")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, numAboveThreshold, [10*i for i in range(1, 10)],
        colors=[(0, i/10.0, 0, 1) for i in range(1, 10)])
plt.contour(xVals, yVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR])

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, fracAboveThreshold, [0.5 * i for i in range(1, 10)],
        colors=[(1, i/10.0, 1-i/10.0, 1) for i in range(1, 10)])
plt.contour(xVals, yVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed')


# Plot thresholds
#plt.plot(L_0_RANGE, [0.62-0.16, 0.62-0.16], c='blue', linewidth=1)
#plt.plot(L_0_RANGE, [0.62+0.15, 0.62+0.15], c='blue', linewidth=1)
#plt.plot([(0.88-0.41) * 1e34, (0.88-0.41) * 1e34], SIGMA_L_RANGE, c='blue', linewidth=1)
#plt.plot([(0.88+0.79) * 1e34, (0.88+0.79) * 1e34], SIGMA_L_RANGE, c='blue', linewidth=1)


plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^')
#plt.scatter(minPoint[0], minPoint[1], c='cyan')


# Observation
shade(numAboveThreshold, NUM_PULSARS_ABOVE_THRESHOLD, xVals, yVals)
shade(fracAboveThreshold, FRAC_ABOVE_THRESHOLD, xVals, yVals, True)


# Final points
if DRAW_PLOEG_POINT:
    plt.plot(ploegPoint[0], ploegPoint[1], markeredgecolor='black', markerfacecolor="C6", marker='s')

custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', linestyle='None'),
                Line2D([], [], markeredgecolor='black', markerfacecolor="C6", marker='s', linestyle='None'),]
plt.legend(custom_lines, ['$N_\\textrm{r} = 47$', '$R_\\textrm{r}=0.2$', "Globular clusters", "GCE"], loc="lower left")
plt.xlim(xVals[0], xVals[-1])
plt.ylim(yVals[0], yVals[-1])

plt.tight_layout()

# Save
if(DRAW_EXTRA_CONTOURS):
    plt.savefig("contour-overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("log-normal-step.pdf")

plt.show()
