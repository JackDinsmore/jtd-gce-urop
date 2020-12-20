"""
The three papers used most so far in this study (Zhong et al. with the power law luminosity function.
Hooper and Linden with the log normal function, and Ploeg et al. with the numerical function), all
use slightly different ranges for the spectrum of the GCE. Zhong et al. uses a range of 0.275 to 
51.9 GeV, whereas both Hooper and Linden and Ploeg et al. use 0.1 to 100 GeV.

Our calculations make use of the total luminosity of the GCE, which is derived from these cutoffs
and the spectrum's distribution. What total GCE luminosities do each of these cutoffs project?

To calculate this, we use a power law distribution for the GCE spectrum with an exponential cutoff
as described by Calore, Cholis, and Weniger 2014 (https://arxiv.org/pdf/1409.0042.pdf). In figure 20,
they display their best fit for the parameters of this distribution to the GCE spectrum. We use a
point in the center of the given contour.
"""

from scipy.special import gamma, gammainc

L_EXCESS = 6.37e36  # All units are in ergs per second
L_CUT = 2.25 # GeV
ALPHA = 0.875
FERMI_LMIN = 0.275
FERMI_LMAX = 51.9

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def integral(a, b): # Integral of the luminosity function
    # Uses units of GeV
    return -(Gamma(1 - ALPHA, b/L_CUT) - Gamma(1 - ALPHA, a / L_CUT))* L_CUT ** (1 - ALPHA) # Integral of the luminosity function
    
def lintegral(a, b): # Integral of the luminosity function times l
    return -(Gamma(2 - ALPHA, b/L_CUT) - Gamma(2 - ALPHA, a / L_CUT))* L_CUT ** (2 - ALPHA)

def getLuminosity(lmin, lmax):
    unscaledFermiLum = lintegral(FERMI_LMIN, FERMI_LMAX)# GeV
    scale = L_EXCESS / unscaledFermiLum# Converts to ergs automatically
    unscaledLum = lintegral(lmin, lmax)
    return unscaledLum * scale

def convertPhotonsPerSecondToErgs(lmin, lmax):
    averageGeVEnergy = lintegral(lmin, lmax) / integral(lmin, lmax)
    return averageGeVEnergy * 0.00160218


fermilab = getLuminosity(FERMI_LMIN, FERMI_LMAX)
lognormal = getLuminosity(0.1, 100)
nptf = getLuminosity(1.893, 11.943)

print("All units are in ergs per second.")
print("GCE luminosity I was using: {0}".format(L_EXCESS))
print()
print("Fermilab GCE luminosity (0.275-51.9 GeV): {0}".format(fermilab))
print("Log normal and Ploeg GCE luminosity (0.1-100 GeV): {0}".format(lognormal))
print("Percent difference: {0}%".format(abs(lognormal - fermilab) / (lognormal + fermilab) * 100))
print()
print("NPTF GCE luminosity (1.893 to 11.943 GeV): {0}".format(nptf))
print()
print("Photon energy in ergs, NPTF: ", convertPhotonsPerSecondToErgs(1.893, 11.943))