from matplotlib import pyplot as plt
from math import pi
from scipy.special import gammainc, gamma
from matplotlib.lines import Line2D
import numpy as np
from numpy import exp, sqrt, log10, log
import ploegload as pl

plt.style.use('latex')

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

CM_PER_KPC_SQUARED = 9.523396e+42
ERGS_PER_PHOTON = 0.00545625167499331
DIST_TO_CENTER = 8.5

ploegFunc = pl.LuminosityFunction(pl.DISK)

NUM_PLOT_POINTS = 200
PLOT_LOG_MIN = 29.0
PLOT_LOG_MAX = ploegFunc.logxs[-1]

ALPHA = 1.94
L_MIN = 1.0e29
L_MAX = 1.0e35
L_THRESH = 34
L0_LOG_NORMAL = 8.8e33
SIGMA_LOG_NORMAL = 0.62
L0_PLOEG_FIT = 1.61e32
SIGMA_PLOEG_FIT = 0.700
N_BELOW_NFW = -0.66
N_ABOVE_NFW = 18.2
L_BREAK_NFW = 1.76e-10 * ERGS_PER_PHOTON * 4 * pi * (DIST_TO_CENTER * DIST_TO_CENTER * CM_PER_KPC_SQUARED)
N_BELOW_DISK = 1.40
N_ABOVE_DISK = 17.5
L_BREAK_DISK = 6.8e-9 * ERGS_PER_PHOTON * 4 * pi * (DIST_TO_CENTER * DIST_TO_CENTER * CM_PER_KPC_SQUARED)


def powerLaw(x):
    return log(10) * x * x**(-ALPHA) * exp(- x / L_MAX) / (Gamma(1 - ALPHA, L_MIN / L_MAX) * L_MAX ** (1 - ALPHA))

def logNormal(x, L0, sigma):
    return log(10) * x * log10(exp(1)) / (sigma * sqrt(2 * pi) * x) * exp(- (log10(x) - log10(L0))**2 / (2 * sigma**2))

def ploeg(x):
    return log(10) * x * ploegFunc(x)

def nptf(x, nBelow, nAbove, LBreak):
    if x < LBreak:
        return log(10) * x * (x / LBreak) ** (-nBelow)
    else:
        return log(10) * x * (x / LBreak) ** (-nAbove)


fig, ax = plt.subplots(figsize=(6, 4))

plt.xlabel("$\log_{10}(L)$")
plt.ylabel("$L\\frac{dN}{dL}$")
plt.title("Luminosity functions")

x = np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, NUM_PLOT_POINTS)

ys = [
    powerLaw(10**x),
    logNormal(10**x, L0_LOG_NORMAL, SIGMA_LOG_NORMAL),
    logNormal(10**x, L0_PLOEG_FIT, SIGMA_PLOEG_FIT),
    [ploeg(10**l) for l in x],
    [nptf(10**l, N_BELOW_NFW, N_ABOVE_NFW, L_BREAK_NFW) for l in x],
    [nptf(10**l, N_BELOW_DISK, N_ABOVE_DISK, L_BREAK_DISK) for l in x],
]
names = ["Power law", "Log normal", "Custom", "Log normal fit", "NFW PS", "Disk PS"]

for i in range(len(ys)):
    plt.plot(x, ys[i] / max(ys[i]), label = names[i])
plt.axvline(x=L_THRESH, color='k')  

plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig("superimposed.png")



fig, ax = plt.subplots(figsize=(6, 4))

plt.xlabel("$\log_{10}(L)$")
plt.ylabel("$L^2\\frac{dN}{dL}$")
plt.title("Luminosity functions")

x = np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, NUM_PLOT_POINTS)

ys = [
    10**x * powerLaw(10**x),
    10**x * logNormal(10**x, L0_LOG_NORMAL, SIGMA_LOG_NORMAL),
    10**x * logNormal(10**x, L0_PLOEG_FIT, SIGMA_PLOEG_FIT),
    [10**l * ploeg(10**l) for l in x],
    [10**l * nptf(10**l, N_BELOW_NFW, N_ABOVE_NFW, L_BREAK_NFW) for l in x],
    [10**l * nptf(10**l, N_BELOW_DISK, N_ABOVE_DISK, L_BREAK_DISK) for l in x],
]
names = ["Power law", "Log normal", "Custom", "Log normal fit", "NFW PS", "Disk PS"]

for i in range(len(ys)):
    plt.plot(x, ys[i] / max(ys[i]), label = names[i])
plt.axvline(x=L_THRESH, color='k')  

plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig("superimposed-times-l.png")

plt.show()