'''
Method:
    By integrating the NFW distribution over space, I'm going to calculate
    the number of particles present in a given ROI and divide by the
    number of particles present in my ROI. I will use this factor to
    convert fluxes, under the assumption that the flux goes like the
    number of sources.
'''

'''Convert a GCE luminosity to a flux, assuming an NFW-squared distributed
pattern of pulsars with delta function luminosity. See Jan summary for
calculation.'''

from math import pi, cos, sqrt

DELTA_ANGLE = 0.001
DELTA_S = 0.01
R_C=8.5
R_S=20
INTEGRAL_LIMIT=100

def NFWSquared(r, gamma):
    return ((r/R_S)**(-gamma) * (1 + r/R_S)**(-3+gamma))**2

def calculateIntegral(maxB, minB, maxL, gamma):# Give angles in units of degrees
    integral = 0
    b = -maxB * pi / 180.0
    while b < maxB * pi / 180.0:
        volumeElement = cos(b) * DELTA_ANGLE * DELTA_ANGLE * DELTA_S
        l = -maxL * pi / 180.0
        if abs(b) < minB * pi / 180.0:
            b += DELTA_ANGLE
            continue
        while l < maxL * pi / 180.0:
            s=0
            while s < INTEGRAL_LIMIT:
                r = sqrt(R_C*R_C + s * s - 2 * R_C * s * cos(b) * cos(l))
                integral += NFWSquared(r, gamma) * volumeElement
                s += DELTA_S
            l += DELTA_ANGLE
        b += DELTA_ANGLE
    return integral

my_ROI = calculateIntegral(20, 2, 20, 1.2)
my_ROI_Abazajian = calculateIntegral(20, 2, 20, 1.0)
square_ROI = calculateIntegral(20, 0, 20, 1.2)
fifteen_ROI = calculateIntegral(15, 0, 15, 1.2)
ten_ROI = calculateIntegral(10, 0, 10, 1.2)
five_ROI = calculateIntegral(5, 0, 5, 1.2)
small_ROI = calculateIntegral(3.5, 0, 3.5, 1.2)
small_ROI_Abazajian = calculateIntegral(3.5, 0, 3.5, 1.0)

print("Gamma:", my_ROI_Abazajian / my_ROI)
print("40x40 ROI factor:", square_ROI / my_ROI)
print("30x30 ROI factor:", fifteen_ROI / my_ROI)
print("20x20 ROI factor:", ten_ROI / my_ROI)
print("10x10 ROI factor:", five_ROI / my_ROI)
print("Abazajian ROI factor, gamma=1.2:", small_ROI / my_ROI)
print("Abazajian ROI factor, gamma=1.0:", small_ROI_Abazajian / my_ROI_Abazajian)
