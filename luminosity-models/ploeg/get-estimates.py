'''
The purpose of this file is to take the Ploeg luminosity function 
and generate predictions of how many pulsars Fermi should have seen,
and how much flux above the threshold it should have seen.
'''

import ploegload
from matplotlib import pyplot as plt

L_EXCESS = 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN = 1e29
L_MAX = 1.0e35

outStr = ""

for i in range(4):
    f = ploegload.LuminosityFunction(i)
    unscaledNumber = f.integrate(minL=L_MIN, maxL=L_MAX)
    unscaledLum = f.lintegrate(minL=L_MIN, maxL=L_MAX)
    unscaledNumberAbove = f.integrate(minL=L_THRESH, maxL=L_MAX)
    unscaledFluxAbove = f.lintegrate(minL=L_THRESH, maxL=L_MAX)
    print(unscaledNumber, unscaledLum, unscaledNumberAbove, unscaledFluxAbove)

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
""".format(f.getName(), totalNumber, totalLum, numberAbove, R)
    print(addStr)
    outStr+=addStr
    f.display()

f = open("get-estimates-output.txt", 'w')
f.write(outStr)
f.close()
plt.xlim(left=L_MIN, right=L_MAX)
plt.axvline(x=L_THRESH, label="Threshold", color='black')
plt.title("Ploeg luminosity functions")
plt.legend()
plt.savefig("get-estimates.png")
plt.show()