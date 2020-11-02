'''
The purpose of this file is to use the Fermilab paper analysis
(arXiv link: https://arxiv.org/pdf/1911.12369.pdf) for varied
values of L_min and L_max, so that we can get the sensitivity
of their algorithm to these values.

The paper gives a reasonable value of alpha_L at [1.2, 1.5], 
which corresponds to about five thousand pulsars. However,
it states that the 4FGL catalog has found enough pulsars
only to admit a value of alpha_L = 1.93, which corresponds
to about 3,000,000 pulsars. These numbers are reproduced here.

Note that as alpha decreases, the number of pulsars also decreases,
at least in the range of alphas discussed above.
'''




from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors

POWER_STEP = 1.11 # 1 is the minimum

ALPHA_L = 1.5#1.93
L_EXCESS = 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN_RANGE=[1.0e28, 1.0e34]
L_MAX_RANGE=[1.0e34, 1.0e36]

dimMin= int(log(L_MIN_RANGE[1]/L_MIN_RANGE[0]) / log(POWER_STEP))
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def getNumPulsars(lMin, lMax):
    lumHard = 1 / (2 - ALPHA_L) * (lMax**(2-ALPHA_L) - lMin**(2-ALPHA_L))
    Ahard = L_EXCESS / lumHard
    Nhard = Ahard / (1 - ALPHA_L) * (lMax**(1-ALPHA_L) - lMin**(1-ALPHA_L))

    lumExp = lMax**(2-ALPHA_L) * Gamma(2-ALPHA_L, lMin / lMax)
    Aexp = L_EXCESS / lumExp
    Nexp = Aexp *  lMax**(1-ALPHA_L) * Gamma(1-ALPHA_L, lMin / lMax)

    return (Nhard, Nexp)


print("Paper values:", getNumPulsars(1e29, 1e35))
print("Delta function values:", getNumPulsars(L_THRESH, L_THRESH*1.0001))


numPulsarsHard = []
numPulsarsExp = []
for i in range(dimMin):
    lineHard = []
    lineExp = []
    for j in range(dimMax):
        nHard, nExp = getNumPulsars(L_MIN_RANGE[0] * POWER_STEP**i, L_MAX_RANGE[0] * POWER_STEP**j)
        lineHard.append(nExp)
        lineExp.append(200 * abs(nExp-nHard)/(nExp+nHard))
    numPulsarsHard.append(lineHard)
    numPulsarsExp.append(lineExp)


lMinVals = [L_MIN_RANGE[0] * POWER_STEP**i for i in range(dimMin)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]

figsize = plt.figaspect(0.38)
fig, (ax1, ax2) = plt.subplots(ncols=2,figsize=figsize)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Lmax")
ax1.set_ylabel("Lmin")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("Lmax")
ax2.set_ylabel("Lmin")
ax1.set_title("Exp cutoff (alpha={0})".format(ALPHA_L))
ax2.set_title("% Difference between hard and exp (alpha={0})".format(ALPHA_L))


c1 = ax1.pcolor(lMaxVals, lMinVals, numPulsarsHard, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsarsHard]),
                   vmax=max([max(v) for v in numPulsarsHard])), cmap='PuBu_r')
plt.colorbar(c1, ax=ax1, extend='max')
c2 = ax2.pcolor(lMaxVals, lMinVals, numPulsarsExp, cmap='PuBu_r')
plt.colorbar(c2, ax=ax2, extend='max')
plt.savefig("tweaks.png")
plt.show()