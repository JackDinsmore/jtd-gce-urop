from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('latex')

POWER_STEP = 1.1 # 1 is the minimum

L_EXCESS = 6.756e36# 6.37e36  # All units are in ergs per second
L_THRESH = 1.0e34
L_MIN = 1e29
ALPHA_RANGE = [1.1, 2.5]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

DRAW_EXTRA_CONTOURS = False
LINE_COLOR = (0.8, 0.3, 0.1)
SHADE_SCALE=25

dimAlpha= 50
dimMax= int(log(L_MAX_RANGE[1]/L_MAX_RANGE[0]) / log(POWER_STEP))

paperPoint = [1e35, 1.94]

def shade(field, threshold, xs, ys, off=False):
    px = []
    py = []
    for x in range(0 if off else 1, SHADE_SCALE, 1):
        inx = int(float(x) / SHADE_SCALE * field.shape[1])
        for y in range(0 if off else 1, SHADE_SCALE, 1):
            iny = int(float(y) / SHADE_SCALE * field.shape[0])
            if field[iny][inx] < threshold:
                fracx = float(x) / SHADE_SCALE * field.shape[1] - inx
                fracy = float(y) / SHADE_SCALE * field.shape[0] - iny
                px.append(xs[inx] + fracx * (xs[inx+1] - xs[inx]))
                py.append(ys[iny] + fracy * (ys[iny+1] - ys[iny]))
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.5)
        

def Gamma(s, x):
    if(s < 0):
        return (Gamma(s+1, x) - x**s * exp(-x))/s
    return gamma(s) * (1-gammainc(s, x))

def getNumPulsars(alpha, lMax):
    lumExp = lMax**(2-alpha) * Gamma(2-alpha, L_MIN / lMax)
    Aexp = L_EXCESS / lumExp
    Nexp = Aexp *  lMax**(1-alpha) * Gamma(1-alpha, L_MIN / lMax)

    return Nexp

def getNumPulsarsAboveThreshold(alpha, lMax):
    lumExp = lMax**(2-alpha) * Gamma(2-alpha, L_MIN / lMax)
    Aexp = L_EXCESS / lumExp
    nAbove = Aexp *  lMax**(1-alpha) * Gamma(1-alpha, L_THRESH / lMax)
    return nAbove

def getFracLumAboveThreshold(alpha, lMax):
    # return (lum above thresh) / totalLum
    fracAbove = Gamma(2-alpha, L_THRESH / lMax) / (Gamma(2-alpha, L_MIN / lMax))
    return fracAbove


print("""Paper point info:
    Coordinates: L0 = {0}, sigma = {1}
    Total num pulsars: {2}
    Num pulsars above threshold: {3}
    Fraction luminosity above threshold: {4}""".format(paperPoint[1], paperPoint[0], getNumPulsars(paperPoint[1], paperPoint[0]), 
    getNumPulsarsAboveThreshold(paperPoint[1], paperPoint[0]), getFracLumAboveThreshold(paperPoint[1], paperPoint[0])))


numPulsars = []
for i in range(dimAlpha):
    line = []
    for j in range(dimMax):
        nExp = getNumPulsars(ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha), L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    numPulsars.append(np.asarray(line))

numAboveThreshold = []
for i in range(dimAlpha):
    line = []
    for j in range(dimMax):
        nExp = getNumPulsarsAboveThreshold(ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha), L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    numAboveThreshold.append(np.asarray(line))

fracAboveThreshold = []
for i in range(dimAlpha):
    line = []
    for j in range(dimMax):
        nExp = getFracLumAboveThreshold(ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha), L_MAX_RANGE[0] * POWER_STEP**j)
        line.append(nExp)
    fracAboveThreshold.append(np.asarray(line))

numPulsars = np.transpose(np.stack(numPulsars))
numAboveThreshold = np.transpose(np.stack(numAboveThreshold))
fracAboveThreshold = np.transpose(np.stack(fracAboveThreshold))



alphaVals = [ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / dimAlpha) for i in range(dimAlpha)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]



fig, ax = plt.subplots(figsize=(6, 4))

plt.yscale("log")
plt.ylabel("$L_{max}$")
plt.xlabel("$\\alpha$")
plt.title("Power law exponential cutoff, step function")

c1 = plt.pcolor(alphaVals, lMaxVals, numPulsars, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in numPulsars]),
                   vmax=max([max(v) for v in numPulsars])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N$")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(alphaVals, lMaxVals, numAboveThreshold, [10*i for i in range(1, 20)], 
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)], linewidths=1)
plt.contour(alphaVals, lMaxVals, numAboveThreshold, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR], linewidths=2, label="$N_r=47$")

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(alphaVals, lMaxVals, fracAboveThreshold, [0.1*i for i in range(1, 15)], 
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)], linewidths=1)
plt.contour(alphaVals, lMaxVals, fracAboveThreshold, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2, label="$R_r=0.2$")

# Observation
shade(numAboveThreshold, NUM_PULSARS_ABOVE_THRESHOLD, alphaVals, lMaxVals)
shade(fracAboveThreshold, FRAC_ABOVE_THRESHOLD, alphaVals, lMaxVals, True)


# Final points 
plt.scatter(paperPoint[1], paperPoint[0], c='purple')

custom_lines = [Line2D([0], [0], color=LINE_COLOR, lw=2),
                Line2D([0], [0], color=LINE_COLOR, lw=2, dashes=(4, 2))]
plt.legend(custom_lines, ["$N_r=47$", "$R_r=0.2$"])
plt.ylim(lMaxVals[0], lMaxVals[-1])
plt.xlim(alphaVals[0], alphaVals[-1])
plt.tight_layout()


if(DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay.png")

plt.show()