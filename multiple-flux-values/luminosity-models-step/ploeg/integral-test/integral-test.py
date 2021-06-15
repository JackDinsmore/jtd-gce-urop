from matplotlib import pyplot as plt
from math import log10, exp, log, sqrt, erf
import numpy as np

plt.style.use('latex')

L_MIN=29.0
L_MAX=35.0
L0 = 1.6084e+32
SIGMA = 0.7003
PI = 3.1415926535
START_BINS = 75
THRESHOLD = 1e34

def getLogNormal(x):
    return log10(exp(1)) / (SIGMA * sqrt(2 * PI) * x) * exp(- (log10(x) - log10(L0))**2 / (2 * SIGMA * SIGMA))

def integrateLogNormal(start):
    return 0.5 * (1 - erf((log10(start)-log10(L0))/(sqrt(2) * SIGMA)))

def lintegrateLogNormal(start):
    return 0.5 * L0 * exp(SIGMA * SIGMA * log(10)*log(10) / 2) * (1 - erf((log10(start) - log10(L0) - SIGMA * SIGMA * log(10)) / (sqrt(2) * SIGMA)))



class LuminosityFunction:
    def __init__(self, n_bins):
        self.xs = 10**np.linspace(29, 35, n_bins)
        self.ys = [getLogNormal(x) for x in self.xs]

    def integrate(self, start):
        i = 0
        while i < len(self.xs) and self.xs[i] < start:
            i+= 1
        assert(i < len(self.xs))
        assert(i > 0)
        frac = (start - self.xs[i - 1]) / (self.xs[i] - self.xs[i - 1])
        startY = self.ys[i - 1] + (self.ys[i] - self.ys[i - 1]) * frac
        integral = 0.5 * (startY + self.ys[i]) * (self.xs[i] - start)
        while i < len(self.xs) - 1:
            integral += 0.5 * (self.ys[i] + self.ys[i + 1]) * (self.xs[i + 1] - self.xs[i])
            i+= 1
        return integral

    def lintegrate(self, start):
        i = 0
        while i < len(self.xs) and self.xs[i] < start:
            i+= 1
        assert(i < len(self.xs))
        assert(i > 0)
        frac = (start - self.xs[i - 1]) / (self.xs[i] - self.xs[i - 1])
        startY = self.ys[i - 1] + (self.ys[i] - self.ys[i - 1]) * frac
        integral = 1 / 6.0 * (-start * start * (self.ys[i] + 2 * startY) + start * self.xs[i] * (startY - self.ys[i]) + self.xs[i] * self.xs[i] * (2 * self.ys[i] + startY))
        while i < len(self.xs) - 1:
            integral += 1 / 6.0 * (self.ys[i] * (-2 * self.xs[i] * self.xs[i] + self.xs[i] * self.xs[i + 1] + self.xs[i + 1] * self.xs[i + 1])
                + self.ys[i + 1] * (-self.xs[i] * self.xs[i] - self.xs[i] * self.xs[i + 1] + 2 * self.xs[i + 1] * self.xs[i + 1]))
            i+= 1   
        return integral


bins = range(20, 150, 2)
starts = 10**np.linspace(L_MIN+0.01, L_MAX-0.01, START_BINS)

integrateData = []
lintegrateData = []

for b in bins:
    integrateLine = []
    lintegrateLine = []
    for s in starts:
        f = LuminosityFunction(b)
        trueIntegrate = integrateLogNormal(s)
        numIntegrate = f.integrate(s)
        integrateLine.append(abs(trueIntegrate - numIntegrate) / trueIntegrate * 100)
        
        trueLintegrate = lintegrateLogNormal(s)
        numLintegrate = f.lintegrate(s)
        lintegrateLine.append(abs(trueLintegrate - numLintegrate) / trueLintegrate * 100)
    integrateData.append(np.asarray(integrateLine))
    lintegrateData.append(np.asarray(lintegrateLine))
integrateData = np.transpose(np.stack(integrateData))
lintegrateData = np.transpose(np.stack(lintegrateData))



coord = (np.argmin(abs(starts - THRESHOLD)), bins.index(100))
print("Percent off for integrate at threshold: {0}. \nPercent off for lintegrate: {1}.".format(integrateData[coord[0]][coord[1]], lintegrateData[coord[0]][coord[1]]))
coord = (0, bins.index(100))
print("Percent off for integrate at min: {0}. \nPercent off for lintegrate: {1}.".format(integrateData[coord[0]][coord[1]], lintegrateData[coord[0]][coord[1]]))



fig, ax1 = plt.subplots(figsize=(7, 5))
c1 = ax1.pcolor(bins, starts, integrateData, cmap = "Greys_r")
ax1.set_title("Number integral deviation")
ax1.set_xlabel("$N$, number of data points")
ax1.set_ylabel("$L_{min}$, starting value")
ax1.set_yscale('log')
ax1.axvline(x=100, label="$N_{bins}$")
ax1.axhline(y=THRESHOLD, label="$L_{th}$")
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("Percent deviation between numerical and analytical result")

plt.tight_layout()
plt.legend()

plt.savefig("num-integral-deviation.png")


fig, ax2 = plt.subplots(figsize=(7, 5))
c2 = ax2.pcolor(bins, starts, lintegrateData, cmap = "Greys_r")
ax2.set_title("Luminosity integral deviation")
ax2.set_xlabel("$N$, number of data points")
ax2.set_ylabel("$L_{min}$, starting value")
ax2.set_yscale('log')
ax2.axvline(x=100, label="$N_{bins}$")
ax2.axhline(y=THRESHOLD, label="$L_{th}$")
cbar = plt.colorbar(c2, extend='max')
cbar.set_label("Percent deviation between numerical and analytical result")

plt.tight_layout()
plt.legend()

plt.savefig("lum-integral-deviation.png")

plt.show()