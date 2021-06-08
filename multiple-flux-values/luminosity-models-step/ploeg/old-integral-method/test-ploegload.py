import ploegload as pl
import numpy as np
from matplotlib import pyplot as plt
from math import exp, erf, log, log10, sqrt

LOG_L_MIN = 30.0
LOG_L_MAX = 34.0
NUM_POINTS = 100

L0 = 1.6084e+32
SIGMA = 0.7003

f = pl.LuminosityFunction(pl.LOG_NORMAL_FIT)
norm = f.integrate()

def logNormLint(logl):
    #return 0.5 * (1 - erf((logl-log10(L0))/(sqrt(2) * SIGMA)))
    return 0.5 * L0 * exp(SIGMA * SIGMA * log(10)*log(10) / 2) * (1 - erf((logl - log10(L0) - SIGMA * SIGMA * log(10)) / (sqrt(2) * SIGMA)))

logls = np.linspace(LOG_L_MIN, LOG_L_MAX, NUM_POINTS)
calcLints = [f.lintegrate(minL=10**logl) / norm for logl in logls]
realLints = [logNormLint(logl) for logl in logls]
print(sum([abs(calcLints[i] - realLints[i]) / realLints[i] for i in range(len(calcLints))]) / len(calcLints))

plt.plot(logls, calcLints, label="My lintegrate funcion")
plt.plot(logls, realLints, label="Log normal fit")
plt.legend()
plt.tight_layout()
plt.savefig("test-lintegrate.png")
plt.show()