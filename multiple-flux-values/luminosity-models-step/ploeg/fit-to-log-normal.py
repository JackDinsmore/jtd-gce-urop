import ploegload
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from math import sqrt, e, pi
import numpy as np

L_MIN = 29
L_MAX = 33

INTEGRATION_LOG_STEP = 0.01

def fitUnfixedNorm(logl, logmean, sigma, norm):
    return norm  / (sigma * sqrt(2 * pi)) * np.exp(-(logl - logmean)**2 / (2 * sigma**2))

def fitLogNormal(functionNumber):

    xdata =[]
    ydata =[]
    f = open("C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-step/ploeg/data/disk.csv")
    lines = f.read().split('\n')
    f.close()
    for line in lines:
        if line == '': continue
        x, y = line.split(",")# log L, dn/dlog L
        xdata.append(float(x))
        ydata.append(float(y))

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
    plt.text(L_MIN, posy-2 * posy/15, "Norm $= {0}$".format(norm))
    plt.legend()
    plt.savefig("fit-to-log-normal.png")
    return logmean, sigma, norm

print(fitLogNormal(ploegload.DISK)) # 32.20641612096041, 0.7003367947758609
plt.show()