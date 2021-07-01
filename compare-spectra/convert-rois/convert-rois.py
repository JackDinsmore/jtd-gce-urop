from math import pi, cos, sqrt
import numpy as np

DELTA_ANGLE = 0.01
DELTA_S = 0.01
R_C=8.5
R_S=20
INTEGRAL_LIMIT=100

def NFWSquared(r, gamma):
    return ((r/R_S)**(-gamma) * (1 + r/R_S)**(-3+gamma))**2

def get_ring(ring, minB, gamma):
    if ring < minB:
        return 0
    integral = 0

    b = 0
    l = ring
    while b < ring:
        volumeElement = cos(b) * DELTA_ANGLE * DELTA_ANGLE * DELTA_S
        if b < minB:
            b += DELTA_ANGLE
            continue
        s=0
        while s < INTEGRAL_LIMIT:
            r = sqrt(R_C*R_C + s * s - 2 * R_C * s * cos(b) * cos(l))
            if r == 0:
                s += DELTA_S
                continue
            integral += NFWSquared(r, gamma) * volumeElement
            s += DELTA_S
        b += DELTA_ANGLE

    b = ring
    l = 0
    volumeElement = cos(b) * DELTA_ANGLE * DELTA_ANGLE * DELTA_S
    while l <= ring:
        s=0
        while s < INTEGRAL_LIMIT:
            r = sqrt(R_C*R_C + s * s - 2 * R_C * s * cos(b) * cos(l))
            if r == 0:
                s += DELTA_S
                continue
            integral += NFWSquared(r, gamma) * volumeElement
            s += DELTA_S
        l += DELTA_ANGLE

    if ring == 0:
        # Do not multiply by four; there is only one zero pixel
        return integral

    # The edges of the quadrant are doubled. Remove them
    volumeElement1 = cos(ring) * DELTA_ANGLE * DELTA_ANGLE * DELTA_S\
        + DELTA_ANGLE * DELTA_ANGLE * DELTA_S
    s=0
    while s < INTEGRAL_LIMIT:
        r = sqrt(R_C*R_C + s * s - 2 * R_C * s * cos(ring))
        if r == 0:
            s += DELTA_S
            continue
        integral -= 0.5 * NFWSquared(r, gamma) * volumeElement
        s += DELTA_S

    return 4 * integral

def get_rings(maxB, minB, gamma):# Give angles in units of degrees
    rings = np.arange(0, maxB * pi / 180, DELTA_ANGLE)
    return rings, [get_ring(r, minB * pi / 180, gamma) for r in rings]

def get_roi(angle, rings, nums):
    integral = 0
    for i, r in enumerate(rings):
        if r < angle * pi / 180:
            integral += nums[i]
    return integral

print("ROI 1/4")
rings_mask, num_mask = get_rings(20, 2, 1.2)
print(num_mask)
print("ROI 2/4")
rings_no_mask, num_no_mask = get_rings(20, 0, 1.2)
print(num_no_mask)
#print("ROI 3/4")
#rings_mask_gamma, num_mask_gamma = get_rings(20, 2, 1.1)
#print("ROI 4/4")
#rings_no_mask_gamma, num_no_mask_gamma = get_rings(5, 0, 1.1)

my_ROI = get_roi(20, rings_mask, num_mask)
print("40x40 ROI factor:", get_roi(20, rings_no_mask, num_no_mask) / my_ROI)
print("30x30 ROI factor:", get_roi(15, rings_no_mask, num_no_mask) / my_ROI)
print("20x20 ROI factor:", get_roi(10, rings_no_mask, num_no_mask) / my_ROI)
print("15x15 ROI factor:", get_roi(7.5, rings_no_mask, num_no_mask) / my_ROI)
print("10x10 ROI factor:", get_roi(5, rings_no_mask, num_no_mask) / my_ROI)
print("7x7 ROI factor:", get_roi(3.5, rings_no_mask, num_no_mask) / my_ROI)

#my_ROI_Abazajian = get_roi(20, rings_mask_gamma, num_mask_gamma)
#print("Gamma:", my_ROI_Abazajian / my_ROI)
#print("7x7 ROI factor, gamma=1.1:", get_roi(3.5, rings_no_mask_gamma, num_no_mask_gamma) / my_ROI_Abazajian)
