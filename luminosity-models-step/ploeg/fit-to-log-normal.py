import ploegload
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from math import sqrt, e, pi
import numpy as np

L_MIN = 30
L_MAX = 35

def fitUnfixedNorm(logl, logmean, sigma, norm):
    return norm  / (sigma * sqrt(2 * pi)) * np.exp(-(logl - logmean)**2 / (2 * sigma**2))

def fitFixedNorm(logl, logmean, sigma):
    return 1  / (sigma * sqrt(2 * pi)) * np.exp(-(logl - logmean)**2 / (2 * sigma**2))

def fitLogNormal(functionNumber, fixNorm):
    f = ploegload.LuminosityFunction(ploegload.DISK)
    xdata = [i[0] for i in f.data if L_MIN < i[0] < L_MAX]
    ydata = [i[1] for i in f.data if L_MIN < i[0] < L_MAX]
    if fixNorm:
        popt, pcov = curve_fit(fitFixedNorm, xdata, ydata, bounds=((L_MIN, 0.1), (L_MAX, 2)))
        logmean, sigma = popt
        fitYData = [fitFixedNorm(l, logmean, sigma) for l in xdata]
    else:
        popt, pcov = curve_fit(fitUnfixedNorm, xdata, ydata, bounds=((L_MIN, 0.1, 0), (L_MAX, 2, 2)))
        logmean, sigma, norm = popt
        fitYData = [fitUnfixedNorm(l, logmean, sigma, norm) for l in xdata]
    posy = max(np.max(ydata), np.max(fitYData))
    plt.figure()
    plt.title("Fitting Ploeg Luminosity function to log normal")
    plt.xlabel("log10(Luminosity)")
    plt.ylabel("Probability density")
    plt.plot(xdata, ydata, label="Original")
    plt.plot(xdata, fitYData, label="Fit")
    plt.text(L_MIN, posy, "$\log_{{10}}(L_0) = {0}$".format(logmean))
    plt.text(L_MIN, posy-posy/15, "$\sigma = {0}$".format(sigma))
    if not fixNorm:
        plt.text(L_MIN, posy-2 * posy/15, "Norm $= {0}$".format(norm))
    plt.legend()
    plt.savefig("fit-to-log-normal.png")
    if fixNorm:
        plt.savefig("fit-to-log-normal-fix-norm.png")
        return logmean, sigma
    else:
        plt.savefig("fit-to-log-normal-unfix-norm.png")
        return logmean, sigma, norm

print(fitLogNormal(ploegload.DISK, True)) # 32.206, 0.70585
print(fitLogNormal(ploegload.DISK, False))
plt.show()