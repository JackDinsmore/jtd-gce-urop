from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import numpy as np
from scipy.special import gammainc, gamma

plt.style.use('latex')

def plotCurves(filename):

    f = open("spectrum-data/" + filename + "-data.csv")
    t = f.readlines()
    f.close

    xdata = []
    ydata = []

    L_MIN = 0.1
    L_MAX = 100.0
    ERG_PER_GEV = 0.00160218
    SIZE_OF_ROI = 0.42882 # steradians
    Y_SCALE = 2.0e8

    for line in t:
        if line == '': continue
        x, y = line.split(', ')
        x, y = (10**float(x), 10**float(y))
        xdata.append(x)
        ydata.append(y * Y_SCALE)
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)

    lbreak, alphaBelow, alphaAbove = 2.06, 1.42, 2.63

    def powerLaw(x, scale):
        alpha = np.full_like(x, alphaAbove)
        alpha[x < lbreak] = alphaBelow
        return x**2 * scale * (x / lbreak)**(-alpha)

    popt, pcov = curve_fit(powerLaw, xdata, ydata, p0=(10), bounds = ((1.0), (100.0)))

    realnorm = popt * (lbreak**alphaBelow / (alphaBelow - 2) * (L_MIN**(2 - alphaBelow) - lbreak**(2 - alphaBelow)) + 
                            lbreak**alphaAbove / (alphaAbove - 2) * (lbreak**(2 - alphaAbove) - L_MAX**(2 - alphaAbove))) / Y_SCALE

    print(filename, "scale:", popt)
    print(filename, "norm (GeV / cm^2 / s / sr):", realnorm)
    print(filename, "norm (erg / cm^2 / s):", realnorm * ERG_PER_GEV * SIZE_OF_ROI)
    fitdata = powerLaw(xdata, popt)

    plt.plot(xdata, ydata / Y_SCALE, label = filename + " data")
    plt.plot(xdata, fitdata / Y_SCALE, label = filename + " fit")

plotCurves("calore")
plotCurves("fermilab")

plt.title("Broken power law fit to GCE flux")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy (GeV)")
plt.ylabel("$E^2 \\frac{dN}{dE}$")
plt.legend()
plt.tight_layout()
plt.savefig("Power law fit.png")
plt.show()