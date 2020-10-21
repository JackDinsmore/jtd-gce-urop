from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors

POWER_STEP = 1.1 #1.05# 1 is the minimum


ALPHA_L_MIN = 1.21
ALPHA_L_MAX = 2.5
PLOT_EVERY = 2
ALPHA_L_INCREMENT = 0.02

ALPHA_L = ALPHA_L_MIN
L_EXCESS = 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN_RANGE=[1.0e28, 1.0e34]
L_MAX_RANGE=[1.0e33, 1.0e37]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/4.0

dimMin= int(log(L_MIN_RANGE[1]/L_MIN_RANGE[0]) / log(POWER_STEP))
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))

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
    fracAbove = Gamma(2-ALPHA_L, L_THRESH / lMax) / (Gamma(2-ALPHA_L, lMin / lMax) - Gamma(2-ALPHA_L, L_THRESH / lMax))
    return fracAbove

fig, ax = plt.subplots(figsize=(6,4))
plt.text(0.95, 0.95, 'Greens: number limit\nReds: luminosity limit', horizontalalignment='right',verticalalignment='top', transform=ax.transAxes, backgroundcolor=(0, 0, 0, 0.3))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Lmax")
plt.ylabel("Lmin")
plt.title("Observational constrants as a function of alpha".format(ALPHA_L))


lMinVals = [L_MIN_RANGE[0] * POWER_STEP**i for i in range(dimMin)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]

intersectionPointsMin = []
intersectionPointsMax = []

plotNum = 0

while ALPHA_L <= ALPHA_L_MAX:
    minErrorCoordinates = (0, 0)
    minError = -1
    numAboveThreshold = []
    for i in range(dimMin):
        line = []
        for j in range(dimMax):
            numAbove = getNumPulsarsAboveThreshold(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
            line.append(numAbove)
        numAboveThreshold.append(line)

    fracAboveThreshold = []
    for i in range(dimMin):
        line = []
        for j in range(dimMax):
            fracAbove = getFracLumAboveThreshold(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
            line.append(fracAbove)
            totalError = ((fracAbove - FRAC_ABOVE_THRESHOLD)/FRAC_ABOVE_THRESHOLD)**2 + ((numAboveThreshold[i][j] - NUM_PULSARS_ABOVE_THRESHOLD)/NUM_PULSARS_ABOVE_THRESHOLD)**2
            if minError > totalError or minError < 0:
                minError = totalError
                minErrorCoordinates = (i, j)
        fracAboveThreshold.append(line)

    if minError < 0.001 and minError > 0:
        # Only consider points which are reasonably close to their goal
        intersectionPointsMin.append(lMinVals[minErrorCoordinates[0]])
        intersectionPointsMax.append(lMaxVals[minErrorCoordinates[1]])

    if plotNum % PLOT_EVERY == 0:
        plt.contour(lMaxVals, lMinVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], 
        colors=[(0, (ALPHA_L-ALPHA_L_MIN)/(ALPHA_L_MAX - ALPHA_L_MIN), 0, 1)], linewidths=1)
        plt.contour(lMaxVals, lMinVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], 
            colors=[(1, (ALPHA_L-ALPHA_L_MIN)/(ALPHA_L_MAX - ALPHA_L_MIN), 1-(ALPHA_L-ALPHA_L_MIN)/(ALPHA_L_MAX - ALPHA_L_MIN), 1)], linewidths=1)

    ALPHA_L += ALPHA_L_INCREMENT
    plotNum+=1

plt.plot(intersectionPointsMax, intersectionPointsMin)
plt.savefig("contour-vary-alpha-exp.png")

plt.show()