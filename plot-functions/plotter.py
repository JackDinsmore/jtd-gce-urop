from matplotlib import pyplot as plt
from math import pi
from scipy.special import gammainc, gamma
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * np.exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

CM_PER_KPC_SQUARED = 9.523396e+42
ERGS_PER_PHOTON = 0.00545625167499331
DIST_TO_CENTER = 8.5

LUM_TO_FLUX = 1.1093417307914119e-46

NUM_PLOT_POINTS = 200
PLOT_LOG_MIN = 29.0
PLOT_LOG_MAX = 35.0

ALPHA = 1.94
L_MIN = 1.0e29
L_MAX = 1.0e35
L_THRESH = 10e34
L0_LOG_NORMAL = 8.8e33
SIGMA_LOG_NORMAL = 0.62
L0_PLOEG_FIT = 1.61e32
SIGMA_PLOEG_FIT = 0.700
L0_GAUTAM_FIT = 3.91983577e+32
SIGMA_GAUTAM_FIT = 0.937184991
N_BELOW_NFW = -0.66
N_ABOVE_NFW = 18.2
L_BREAK_NFW = 1.76e-10 * ERGS_PER_PHOTON * 4 * pi * (DIST_TO_CENTER * DIST_TO_CENTER * CM_PER_KPC_SQUARED)
N_BELOW_DISK = 1.40
N_ABOVE_DISK = 17.5
L_BREAK_DISK = 6.8e-9 * ERGS_PER_PHOTON * 4 * pi * (DIST_TO_CENTER * DIST_TO_CENTER * CM_PER_KPC_SQUARED)


def powerLaw(x, lmin, lmax, alpha):
    return x**-alpha * np.exp(-x / lmax) / (Gamma(1-alpha, lmin/lmax) * lmax**(1-alpha))
def logNormal(x, L0, sigma):
    return np.log10(np.exp(1)) / (sigma * np.sqrt(2 * np.pi) * x)* np.exp(-(np.log10(x) - np.log10(L0))**2 / (2 * sigma**2))
def nptf(x, n1, n2, LBreak):
    if x < LBreak:
        return (1 - n1) * (1 - n2) / (LBreak * (n1 - n2)) * (x/LBreak)**(-n1)
    else:
        return (1 - n1) * (1 - n2) / (LBreak * (n1 - n2)) * (x/LBreak)**(-n2)
def lum_to_flux(lum):
    return lum * LUM_TO_FLUX
def flux_to_lum(flux):
    return flux / LUM_TO_FLUX

fig, ax1=plt.subplots()

ax1.set_xlabel("Luminosity $L$ [erg / s]")
ax1.set_ylabel("$\\frac{dN}{dL}$")
ax1.set_xscale('log')

x = 10**np.linspace(PLOT_LOG_MIN, PLOT_LOG_MAX, NUM_PLOT_POINTS)

ys = [
    powerLaw(x, L_MIN, L_MAX, ALPHA),
    logNormal(x, L0_LOG_NORMAL, SIGMA_LOG_NORMAL),
    logNormal(x, L0_PLOEG_FIT, SIGMA_PLOEG_FIT),
    logNormal(x, L0_GAUTAM_FIT, SIGMA_GAUTAM_FIT),
    [nptf(l, N_BELOW_NFW, N_ABOVE_NFW, L_BREAK_NFW) for l in x],
    [nptf(l, N_BELOW_DISK, N_ABOVE_DISK, L_BREAK_DISK) for l in x],
]
names = ["Power law", "Log normal GCL", "Log normal GCE", "Log normal AIC", "NPTF NFW PS", "NPTF Disk PS"]
styles = ["solid", "dashed", "dashed", "dashed", "dotted", "dotted"]
colors = ["blue", "orange", "purple", "red", "fuchsia", "slategray"]

for i in range(len(ys)):
    ax1.plot(x, np.abs(ys[i]) / np.max(np.abs(ys[i])), label = names[i], linestyle=styles[i],
        c=colors[i])
#ax1.axvline(x=L_THRESH, color='k')

ax2 = ax1.secondary_xaxis('top', functions=(lum_to_flux, flux_to_lum))
ax2.set_xlabel("Flux $F$ [erg / cm$^2$ / s]")

ax1.legend(loc='upper left')
plt.tight_layout()

plt.savefig("lum-funcs.pdf")
plt.show()



fig, ax1 = plt.subplots()

ax1.set_xlabel("Luminosity $L$ [erg / s]")
ax1.set_ylabel("$L\\frac{dN}{dL}$")
ax1.set_xscale('log')

ys = [x * y for y in ys]

for i in range(len(ys)):
    ax1.plot(x, np.abs(ys[i]) / np.max(np.abs(ys[i])), label = names[i], linestyle=styles[i],
    c=colors[i])
#ax1.axvline(x=L_THRESH, color='k')

ax2 = ax1.secondary_xaxis('top', functions=(lum_to_flux, flux_to_lum))
ax2.set_xlabel("Flux $F$ [erg / cm$^2$ / s]")

ax1.legend(loc='upper left')
plt.tight_layout()

plt.savefig("l-lum-funcs.pdf")
plt.show()





fig, ax1 = plt.subplots()

ax1.set_xlabel("Luminosity $L$ [erg / s]")
ax1.set_ylabel("$L^2\\frac{dN}{dL}$")
ax1.set_xscale('log')

ys = [x * y for y in ys]

for i in range(len(ys)):
    ax1.plot(x, np.abs(ys[i]) / np.max(np.abs(ys[i])), label = names[i], linestyle=styles[i],
    c=colors[i])
#ax1.axvline(x=L_THRESH, color='k')
ax2 = ax1.secondary_xaxis('top', functions=(lum_to_flux, flux_to_lum))
ax2.set_xlabel("Flux $F$ [erg / cm$^2$ / s]")

ax1.legend(loc='upper left')
plt.tight_layout()

plt.savefig("l2-lum-funcs.pdf")
plt.show()
