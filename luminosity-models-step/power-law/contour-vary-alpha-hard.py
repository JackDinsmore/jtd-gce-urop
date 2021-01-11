from matplotlib import pyplot as plt
from math import log, exp
import matplotlib.colors as colors

POWER_STEP = 1.1#1.05 # 1 is the minimum


ALPHA_L_MIN = 1.21
ALPHA_L_MAX = 2.5
PLOT_EVERY = 2#4
ALPHA_L_INCREMENT = 0.02

ALPHA_L = ALPHA_L_MIN
L_EXCESS = 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN_RANGE=[1.0e28, 1.0e34]
L_MAX_RANGE=[1.0e34, 1.0e36]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

dimMin= int(log(L_MIN_RANGE[1]/L_MIN_RANGE[0]) / log(POWER_STEP))
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))
    
def getNumPulsars(lMin, lMax):
    lum = 1/(ALPHA_L - 2) * (lMin **(2 - ALPHA_L) - lMax**(2 - ALPHA_L))
    A = L_EXCESS / lum
    N = A *  1/(ALPHA_L - 1) * (lMin **(1 - ALPHA_L) - lMax**(1 - ALPHA_L))
    return N

def getNumPulsarsAboveThreshold(lMin, lMax):
    lum = 1/(ALPHA_L - 2) * (lMin **(2 - ALPHA_L) - lMax**(2 - ALPHA_L))
    A = L_EXCESS / lum
    nAbove = A *  1/(ALPHA_L - 1) * (L_THRESH **(1 - ALPHA_L) - lMax**(1 - ALPHA_L))
    return nAbove

def getFracLumAboveThreshold(lMin, lMax):
    # return (lum above thresh) / (total lum)
    fracAbove = (1/(ALPHA_L - 2) * (L_THRESH **(2 - ALPHA_L) - lMax**(2 - ALPHA_L))) / (1/(ALPHA_L - 2) * (lMin **(2 - ALPHA_L) - lMax**(2 - ALPHA_L)))
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
plt.savefig("contour-vary-alpha-hard.png")

plt.show()