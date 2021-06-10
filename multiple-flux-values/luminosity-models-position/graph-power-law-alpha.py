from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('revtex')


FERMILAB_GNFW_FLUX = 5.716826056794996e-10
CALORE_FLUX = 1.3592975571782718e-09
DI_MAURO_FLUX = 1.2953417255755896e-09
AJELLO_FLUX = 1.8268259089431682e-09
FERMILAB_NFW_FLUX = 7.770902203635269e-10
ABAZAJIAN_FLUX = 4.887121798902831e-10

PLOT_SIZE = 50

L_MIN = 1e29
ALPHA_RANGE = [1.1, 2.5]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1e35, 1e29]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.2
TOTAL_FLUX =1.2953417255755896e-09

DRAW_EXTRA_CONTOURS = True
SHOW_NUMBERS = False
LINE_COLOR = (0.8, 0.3, 0.1)
PATH_TO_FILE = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/data/power-law-alpha/"
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
    plt.scatter(px, py, marker=('|' if off else '_'), c='k', sizes = (20,), alpha=0.5)

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

dimMin = len(totalNum)
dimMax = len(totalNum[0])

alphaVals = [ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / PLOT_SIZE) for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots(figsize=(6, 4))

#plt.xscale("log")
plt.yscale("log")
plt.xlabel("$\\alpha$")
plt.ylabel("$L_\mathrm{max}$")
plt.title("Power-law, position-dependent")

if SHOW_NUMBERS:
    c1 = plt.pcolor(alphaVals, lMaxVals, totalNum, 
                    norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                    vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
    cbar = plt.colorbar(c1, extend='max')
    cbar.set_label("$N$")

# Greens
if DRAW_EXTRA_CONTOURS:
    plt.contour(alphaVals, lMaxVals, numSeen * FERMILAB_GNFW_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["blue"], linewidths=1, label="$N_r$ Fermilab gNFW")
    plt.contour(alphaVals, lMaxVals, numSeen * ABAZAJIAN_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["green"], linewidths=1, label="$N_r$ Abazajian")
    plt.contour(alphaVals, lMaxVals, numSeen * AJELLO_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=["teal"], linewidths=1, label="$N_r$ Ajello")
plt.contour(alphaVals, lMaxVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR], linewidths=2, label="$N_r=47$")

# Reds
plt.contour(alphaVals, lMaxVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2, label="$R_r=0.2$")

# Observation
if DRAW_EXTRA_CONTOURS:
    shade(numSeen * AJELLO_FLUX / DI_MAURO_FLUX, NUM_PULSARS_ABOVE_THRESHOLD, alphaVals, lMaxVals)
else:
    shade(numSeen, NUM_PULSARS_ABOVE_THRESHOLD, alphaVals, lMaxVals)
shade(lumSeen, FRAC_ABOVE_THRESHOLD, alphaVals, lMaxVals, True)


# Final points

plt.scatter(paperPoint[0], paperPoint[1], c='purple')


custom_lines = [Line2D([0], [0], color=LINE_COLOR, lw=2),
                Line2D([0], [0], color="blue", lw=1),
                Line2D([0], [0], color="green", lw=1),
                Line2D([0], [0], color="teal", lw=1),
                Line2D([0], [0], color=LINE_COLOR, lw=2, dashes=(4, 2)),]
plt.legend(custom_lines, ['$N_r$ Di Mauro', "$N_r$ Fermilab gNFW", "$N_r$ Abazajian", "$N_r$ Ajello", '$R_r$'], loc="upper right")
plt.ylim(lMaxVals[0], lMaxVals[-1])
plt.xlim(alphaVals[0], alphaVals[-1])
plt.tight_layout()


if(DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay-alpha-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("overlay-alpha.png")

plt.show()