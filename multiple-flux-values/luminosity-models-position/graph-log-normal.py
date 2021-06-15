from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
import numpy as np

from matplotlib.lines import Line2D

plt.style.use('revtex-presentation')


PLOT_SIZE = 50
L_0_RANGE=[1.0e30, 2.0e36]#[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]

MULTIPLIER = 1

lOPowerStep = (L_0_RANGE[1] / L_0_RANGE[0]) ** (1.0 / PLOT_SIZE)

FERMILAB_GNFW_FLUX = 5.716826056794996e-10
CALORE_FLUX = 1.3592975571782718e-09
DI_MAURO_FLUX = 1.2953417255755896e-09
AJELLO_FLUX = 1.8268259089431682e-09
FERMILAB_NFW_FLUX = 7.770902203635269e-10
ABAZAJIAN_FLUX = 4.887121798902831e-10

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0
TOTAL_FLUX = 7.494712733226778e-10

DRAW_EXTRA_CONTOURS = True
LINE_COLOR = (0.8, 0.3, 0.1)
DRAW_PLOEG_POINT = True

paperPoint = [0.88e34, 0.62]
ploegPoint = [10**32.206, 0.70585]
SHOW_NUMBERS = False

PATH_TO_FILE = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/data/log-normal/"
SHADE_SCALE=25

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
    plt.scatter(px, py, marker=('|' if off else '_'), c='black', sizes = (20,), alpha=0.5)

# ========================== Load data ===========================

totalNum = []
numSeen = []
lumSeen = []

f = open(PATH_TO_FILE + "total-num.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    totalNum.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "num-seen.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    numSeen.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "lum-seen.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item) / TOTAL_FLUX)
    lumSeen.append(np.asarray(enterLine))

totalNum = np.transpose(np.stack(totalNum, axis=0))
numSeen = np.transpose(np.stack(numSeen, axis=0))
lumSeen = np.transpose(np.stack(lumSeen, axis=0))

# ========================= Display data =========================

xVals = [L_0_RANGE[0] * lOPowerStep**i for i in range(PLOT_SIZE)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / PLOT_SIZE * j for j in range(PLOT_SIZE)]


fig, ax = plt.subplots()
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])

plt.xscale("log")
plt.ylabel("$\sigma$")
plt.xlabel("$L_0$")
plt.title("Log-normal, position-dependent")

if SHOW_NUMBERS:
    cols = colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                    vmax=max([max(v) for v in totalNum]))
    c1 = plt.pcolor(xVals, yVals, totalNum, 
                    norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                    vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
    cbar = plt.colorbar(c1, extend='max')
    cbar.set_label("$N$")


# Greens
if DRAW_EXTRA_CONTOURS:
    plt.contour(xVals, yVals, numSeen * FERMILAB_GNFW_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["blue"], linewidths=1, label="$N_r$ Fermilab gNFW")
    plt.contour(xVals, yVals, numSeen * ABAZAJIAN_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["green"], linewidths=1, label="$N_r$ Abazajian")
    plt.contour(xVals, yVals, numSeen * AJELLO_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["teal"], linewidths=1, label="$N_r$ Ajello")
plt.contour(xVals, yVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR], linewidths=2, label="$N_r$")

# Reds
plt.contour(xVals, yVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2, label="$R_r$")

plt.scatter(paperPoint[0], paperPoint[1], c='blue')
#plt.scatter(minPoint[0], minPoint[1], c='cyan')

# Observation
if DRAW_EXTRA_CONTOURS:
    shade(numSeen * AJELLO_FLUX / DI_MAURO_FLUX, NUM_PULSARS_ABOVE_THRESHOLD, xVals, yVals)
else:
    shade(numSeen, NUM_PULSARS_ABOVE_THRESHOLD, xVals, yVals)
shade(lumSeen, FRAC_ABOVE_THRESHOLD, xVals, yVals, True)


# Final points 

if DRAW_PLOEG_POINT:
    plt.scatter(ploegPoint[0], ploegPoint[1], c='green')

custom_lines = [Line2D([0], [0], color=LINE_COLOR, lw=2),
                Line2D([0], [0], color="blue", lw=1),
                Line2D([0], [0], color="green", lw=1),
                Line2D([0], [0], color="teal", lw=1),
                Line2D([0], [0], color=LINE_COLOR, lw=2, dashes=(4, 2)),]
plt.legend(custom_lines, ['$N_r$ Di Mauro', "$N_r$ Fermilab gNFW", "$N_r$ Abazajian", "$N_r$ Ajello", '$R_r$'], loc="lower left")
plt.xlim(xVals[0], xVals[-1])
plt.ylim(yVals[0], yVals[-1])

plt.tight_layout()

# Save
plt.savefig("overlay-log-normal.png")

plt.show()