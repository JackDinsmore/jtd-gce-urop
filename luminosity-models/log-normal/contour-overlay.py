from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors

DIM_TRIALS=100

L_EXCESS = 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_0_RANGE=[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]
powerStep =(L_0_RANGE[1] / L_0_RANGE[0])**(1/DIM_TRIALS)

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

DRAW_EXTRA_CONTOURS = True

def getNumPulsars(L, sigma):
    totalLum = exp(0.5 * sigma**2 * log(10)**2) * L
    A = L_EXCESS / totalLum
    return A * 1 # A times number of pulsars per unit

def getNumPulsarsAboveThreshold(L, sigma):
    totalLum = exp(0.5 * sigma**2 * log(10)**2) * L
    A = L_EXCESS / totalLum
    return A * 0.5 * erfc((log(L_THRESH) - log(L)) / (sqrt(2) * sigma * log(10)))# A times number of pulsars above threshold per unit

def getFracLumAboveThreshold(L, sigma):
    # return (lum above thresh) / (total lum)
    erfArgument = (sqrt(2) * (sigma**2 * log(10)**2 - log(L_THRESH) + log(L))) / (sigma * log(100))
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

def getPaperPointInfo(L0 = 0.88e34, sigma = 0.62):
    print("""Paper point info:
    Coordinates: L0 = {0}, sigma = {1}
    Total num pulsars: {2}
    Num pulsars above threshold: {3}
    Fraction luminosity above threshold: {4}""".format(L0, sigma, getNumPulsars(L0, sigma), 
    getNumPulsarsAboveThreshold(L0, sigma), getFracLumAboveThreshold(L0, sigma)))


getMinNumPulsarsInTriangle()
minPoint = getMinPulsarsWithinOneStdevOfSigma()
paperPoint = [0.88e34, 0.62]
getPaperPointInfo()


numPulsars = []
for j in range(DIM_TRIALS):
    line = []
    for i in range(DIM_TRIALS):
        n = getNumPulsars(L_0_RANGE[0] * powerStep**i, 
            SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j)
        line.append(n)
    numPulsars.append(line)

numAboveThreshold = []
for j in range(DIM_TRIALS):
    line = []
    for i in range(DIM_TRIALS):
        n = getNumPulsarsAboveThreshold(L_0_RANGE[0] * powerStep**i, 
            SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j)
        line.append(n)
    numAboveThreshold.append(line)

fracAboveThreshold = []
for j in range(DIM_TRIALS):
    line = []
    for i in range(DIM_TRIALS):
        frac = getFracLumAboveThreshold(L_0_RANGE[0] * powerStep**i, 
            SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j)
        line.append(frac)
    fracAboveThreshold.append(line)


xVals = [L_0_RANGE[0] * powerStep**i for i in range(DIM_TRIALS)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j for j in range(DIM_TRIALS)]


fig, ax = plt.subplots(figsize=(6, 4))
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])
plt.text(0.95, 0.95, 'Greens: number limit\nReds: luminosity limit\nBold: Fermilab observations', 
    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, color='white', backgroundcolor=(0, 0, 0, 0.3))

plt.xscale("log")
plt.ylabel("Sigma")
plt.xlabel("L_0")
plt.title("Log normal luminosity function")

c1 = plt.pcolor(xVals, yVals, numPulsars, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='PuBu_r')
plt.colorbar(c1, extend='max')

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, numAboveThreshold, [10*i for i in range(1, 10)], 
        colors=[(0, i/10.0, 0, 1) for i in range(1, 10)], linewidths=1)
plt.contour(xVals, yVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[(0, 1, 0)], linewidths=2)

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, fracAboveThreshold, [0.5 * i for i in range(1, 10)], 
        colors=[(1, i/10.0, 1-i/10.0, 1) for i in range(1, 10)], linewidths=1)
plt.contour(xVals, yVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], colors=[(1, 0, 0)], linewidths=2)


# Plot thresholds
plt.plot(L_0_RANGE, [0.62-0.16, 0.62-0.16], c='purple')
plt.plot(L_0_RANGE, [0.62+0.15, 0.62+0.15], c='purple')
plt.plot([(0.88-0.41) * 1e34, (0.88-0.41) * 1e34], SIGMA_L_RANGE, c='purple')
plt.plot([(0.88+0.79) * 1e34, (0.88+0.79) * 1e34], SIGMA_L_RANGE, c='purple')

plt.scatter(paperPoint[0], paperPoint[1], c='purple')
plt.scatter(minPoint[0], minPoint[1], c='black')

# Save
if(DRAW_EXTRA_CONTOURS):
    plt.savefig("contour-overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("contour-overlay.png")


plt.show()