from math import log, exp
from scipy.special import gammainc, gamma
import numpy as np

ALPHA_L = 1.94
L_EXCESS = 1.2953417255755896e-09 / 8.331593765023139e-47#6.756e36# 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN = 1e29
L_MAX = 1e35
MIN_FLUX = 1e-12
MAX_FLUX = 1e-9

FLUX_TO_LUM_RATIO = 1.1093417307914119e-46

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

def fluxToLum(flux):
    return flux / FLUX_TO_LUM_RATIO

def lumToFlux(lum):
    return lum * FLUX_TO_LUM_RATIO

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def getNumPulsars(lMin, lMax, lum):
    lumExp = lMax**(2-ALPHA_L) * Gamma(2-ALPHA_L, lMin / lMax)
    Aexp = L_EXCESS / lumExp

    if lum > L_THRESH:
        return Aexp * lum**-ALPHA_L * np.exp(-lum/L_MAX)
    else:
        return 0

output = open("power-law-step.txt", 'w')
binEdges = 10**np.linspace(start=np.log10(fluxToLum(MIN_FLUX)), stop=np.log10(fluxToLum(MAX_FLUX)), num=41)
bins = [(binEdges[i+1] + binEdges[i])/2 for i in range(len(binEdges) - 1)]
binWidths = [binEdges[i+1] - binEdges[i] for i in range(len(binEdges) - 1)]
sumCounts = 0

for lum in bins:
    output.write(str(lumToFlux(lum)) + ", ")
output.write('\n')
for (i, lum) in enumerate(bins):
    thisCount = getNumPulsars(L_MIN, L_MAX, lum) * binWidths[i]
    sumCounts += thisCount
    output.write(str(thisCount) + ', ')

print(sumCounts)

output.close()