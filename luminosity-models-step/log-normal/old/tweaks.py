'''
The purpose of this file is to use the log-normal luminosity function
(link: https://www.osti.gov/servlets/purl/1305131) for varied
values of sigma and L_0, so that we can get the sensitivity
of the luminosity function to its parameters. It is slightly more realistic
than the power law, because it has a peak. That is why we are trying it.

'''




from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors

DIM_TRIALS=100

L_EXCESS = 6.37e36  # All units are in ergs per second
L_0_RANGE=[1.0e32, 1.0e36]
SIGMA_L_RANGE=[0.001, 3]
powerStep =(L_0_RANGE[1] / L_0_RANGE[0])**(1/DIM_TRIALS)

def getNumPulsars(L, sigma):
    totalLum = exp(0.5 * sigma**2 * log(10)**2) * L
    A = L_EXCESS / totalLum
    return A * 1


print("Paper values:", getNumPulsars(0.88e34, 0.62))


numPulsars = []
for j in range(DIM_TRIALS):
    line = []
    for i in range(DIM_TRIALS):
        n = getNumPulsars(L_0_RANGE[0] * powerStep**i, 
        SIGMA_L_RANGE[0]  + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j)

        line.append(n)
    numPulsars.append(line)


xVals = [L_0_RANGE[0] * powerStep**i for i in range(DIM_TRIALS)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_TRIALS * j for j in range(DIM_TRIALS)]

plt.figure(figsize=(10,7))

plt.xscale("log")
plt.ylabel("Sigma")
plt.xlabel("L_0")
plt.title("Number of MSPs, log normal lum. func.")


c1 = plt.pcolor(xVals, yVals, numPulsars, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='PuBu_r')
plt.colorbar(c1, extend='max')
plt.savefig("tweaks.png")
plt.show()