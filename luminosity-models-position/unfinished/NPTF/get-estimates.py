'''
The purpose of this file is to take the NPTF luminosity function 
and generate predictions of how many pulsars Fermi should have seen,
and how much flux above the threshold it should have seen. We expect
the luminosity to have a short and strong peak below the threshold,
leading to very few visible pulsars and small visible luminosity. 
The paper is found here https://arxiv.org/pdf/1506.05124.pdf,
with a luminosity function described in table 1.
'''

from matplotlib import pyplot as plt
import numpy as np
from math import log10

plt.style.use('latex')

L_EXCESS = 3.3440751340173124e+36  # 1.893 to 11.943 GeV
DIST_FROM_GALACTIC_CENTER_KILOPARSEC = 8.5
CM_PER_KILOPARSEC = 3.086e21
L_THRESH = 1.0e34
ERGS_PER_PHOTON = 0.00545625167499331
PI = 3.1415926535

L_MIN = None
L_MAX = None

global outStr
outStr = ""

class LuminosityFunction:
    def __init__(self, name, nBelow, nAbove, fluxBreak):
        self.name = name
        self.nBelow = nBelow
        self.nAbove = nAbove
        self.lBreak = fluxBreak * ERGS_PER_PHOTON * (4 * PI * (DIST_FROM_GALACTIC_CENTER_KILOPARSEC * CM_PER_KILOPARSEC)**2)
    
    def integrate(self, minL, maxL):
        lowC = self.lBreak if maxL == None else min(maxL, self.lBreak)
        highC = self.lBreak if minL == None else max(minL, self.lBreak)
        low=0
        high=0
        if minL == None:
            #assert(self.nBelow < 1)
            low = self.lBreak / (1 - self.nBelow) * (lowC / self.lBreak)**(1 - self.nBelow)
        elif minL < self.lBreak:
            low = self.lBreak / (1 - self.nBelow) * ((lowC / self.lBreak)**(1 - self.nBelow) - (minL / self.lBreak)**(1 - self.nBelow))
        if maxL == None:
            #assert(self.nAbove > 1)
            high = self.lBreak / (self.nAbove - 1) * (self.lBreak / highC)**(self.nAbove - 1)
        elif maxL > self.lBreak:
            high = self.lBreak / (self.nAbove - 1) * ((self.lBreak / highC)**(self.nAbove - 1) - (self.lBreak / maxL)**(self.nAbove - 1))
        return low+high

    def lintegrate(self, minL, maxL):
        lowC = self.lBreak if maxL == None else min(maxL, self.lBreak)
        highC = self.lBreak if minL == None else max(minL, self.lBreak)
        low=0
        high=0
        if minL == None:
            #assert(self.nBelow < 2)
            low = self.lBreak**2 / (2 - self.nBelow) * (lowC / self.lBreak)**(2 - self.nBelow)
        elif minL < self.lBreak:
            low = self.lBreak**2 / (2 - self.nBelow) * ((lowC / self.lBreak)**(2 - self.nBelow) - (minL / self.lBreak)**(2 - self.nBelow))
        if maxL == None:
            #assert(self.nAbove > 2)
            high = self.lBreak**2 / (self.nAbove - 2) * (self.lBreak / highC)**(self.nAbove - 2)
        elif maxL > self.lBreak:
            high = self.lBreak**2 / (self.nAbove - 2) * ((self.lBreak / highC)**(self.nAbove - 2) - (self.lBreak / maxL)**(self.nAbove - 2))
        return low+high

    def getValue(self, l):
        scale = 1
        if self.nBelow > 0:
            lMin = 1e29 if (L_MIN==None or L_MIN < 1e29) else L_MIN
            scale = (lMin / self.lBreak)**(self.nBelow) # Make the highest point in the luminosity function have value one. Affects only the plot.
        if l < self.lBreak:
            return scale * (l / self.lBreak)**(-self.nBelow)
        else:
            return scale * (l / self.lBreak)**(-self.nAbove)

    def printEstimate(self):
        global outStr
        unscaledNumber = self.integrate(minL=L_MIN, maxL=L_MAX)
        unscaledLum = self.lintegrate(minL=L_MIN, maxL=L_MAX)
        unscaledNumberAbove = self.integrate(minL=L_THRESH, maxL=L_MAX)
        unscaledFluxAbove = self.lintegrate(minL=L_THRESH, maxL=L_MAX)

        scale = L_EXCESS / unscaledLum

        totalNumber = unscaledNumber * scale
        totalLum = unscaledLum * scale # Should be L_EXCESS
        numberAbove = unscaledNumberAbove * scale
        R = unscaledFluxAbove / unscaledLum

        addStr = """{0} luminosity function:
    Total number of pulsars:\t\t{1}
    Total luminosity:\t\t\t{2}
    Number of pulsars above threshold:\t{3}\t(Fermi: 47)
    Fraction of lum above threshold:\t{4}\t(Fermi: 0.2)
""".format(self.name, totalNumber, totalLum, numberAbove, R)
        print(addStr)
        outStr+=addStr

    def display(self):
        x=10**np.linspace(start=(29 if (L_MIN==None or L_MIN < 1e29) else log10(L_MIN)),
            stop=(35 if (L_MAX == None or L_MAX > 1e35) else log10(L_MAX)), num=100)
        y=[self.getValue(l) for l in x]
        plt.plot(x, y, label=self.name)

lumFuncs = [LuminosityFunction("NFW PS", -0.66, 18.2, 1.76e-10),
            LuminosityFunction("Disk PS", 1.40, 17.5, 6.8e-9)
            ]

for l in lumFuncs:
    l.printEstimate()
    l.display()

f = open("get-estimates-output.txt", 'w')
f.write(outStr)
f.close()

plt.axvline(x=L_THRESH, label="Threshold", color='black')
plt.xscale('log')
#plt.yscale('log')
plt.title("NPTF luminosity function")
plt.savefig("get-estimates.png")
plt.xlabel("Luminosity (ergs/s)")
plt.ylabel("$\\frac{dN}{dL}$ (unnormalized)")
plt.legend()

plt.savefig("luminosity-func.png")
plt.show()